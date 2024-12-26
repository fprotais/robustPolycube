#include "blockExtraction.h"
#include <algorithm>
#include <utils/trace.h>
#include "tools.h"

#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;
using namespace rb_data_structure;


void compute_charts(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, const UM::CellFacetAttribute<int>& cfflag, rb_data_structure::Sorted_charts& charts) {
	Triangles surf;
	std::vector<int> surf2cf;
	std::vector<int> cf_is_bnd(m.ncells() * 4);

	surf.points.data->assign(m.points.begin(), m.points.end());
	FOR(cf, 4 * m.ncells()) if (cfflag[cf] != -1) {
		cf_is_bnd[cf] = true;
		int f = surf.create_facets(1);
		surf2cf.push_back(cf);
		FOR(fv, 3) surf.vert(f, fv) = m.facet_vert(cf / 4, cf % 4, fv);
	}

	DisjointSet ds(3 * m.nverts());
	FOR(f, surf.nfacets()) FOR(fv, 3) {
		int dim = (cfflag[surf2cf[f]] / 2);
		ds.merge(dim * surf.nverts() + surf.vert(f, fv), dim * surf.nverts() + surf.vert(f, (fv + 1) % 3));
	}

	std::vector<int> ids;
	int nb_ids = ds.get_sets_id(ids);
	std::vector<int> bnd_id(nb_ids, -1);

	FOR(cf, 4 * m.ncells()) if (cf_is_bnd[cf]) {
		int dim = (cfflag[cf] / 2);
		int rep_v = m.facet_vert(cf/4, cf%4, 0);
		int id = ids[dim * m.nverts() + rep_v];
		if (bnd_id[id] == -1) {
			bnd_id[id] = charts.size();
			charts.charts.push_back(Chart());
			charts.charts.back().dim = dim;
			charts.charts.back().value = polycuboid.points[rep_v][dim];
			charts.charts.back().tetFacets.push_back(cf);
		}
		else {
			charts[bnd_id[id]].tetFacets.push_back(cf);
		}
	}

	std::sort(charts.charts.begin(), charts.charts.end(), [](const Chart& a, const Chart& b) {
		if (a.dim == b.dim) {
			return a.value < b.value;
		}
		else return a.dim < b.dim;
		});

	charts.start_of_dim = { (int)charts.size(), (int)charts.size(), (int)charts.size() };
	FOR(c, charts.size()) if (charts.start_of_dim[charts[c].dim] > c) charts.start_of_dim[charts[c].dim] = c;
	
	charts.cf2chart.resize(m.ncells() * 4, -1);
	FOR(c, charts.size()) for (int cf : charts[c].tetFacets) charts.cf2chart[cf] = c;

	FOR(c, charts.size()) charts[c].id = c;

}

inline bool check_mesh_is_connexe(Hexahedra& m) {
	DisjointSet ds(m.ncells());
	OppositeFacet adj(m);
	FOR(c, m.ncells()) FOR(cf, 6) if (adj.adjacent[6*c+cf] != -1) {
		ds.merge(c, adj.adjacent[6 * c + cf] / 6);
	}
	std::vector<int> tmp;
	int nb_connected_set = ds.get_sets_id(tmp);
	return nb_connected_set == 1;
}

inline bool check_bnd_mesh_is_variety(Hexahedra& m) {
	Trace::warning("Pouet, function check_bnd_mesh_is_variety need implementation.");
	return true;
}

void chart_voxelisation(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, const rb_data_structure::Sorted_charts& charts, rb_data_structure::Block_decomposition& blocks) {

	std::vector<std::array<int, 3>> v_on_chart(m.nverts(), { { -1,-1,-1} });
	for (Chart c : charts.charts) for (int cf : c.tetFacets) FOR(cfv, 3) v_on_chart[m.facet_vert(cf / 4, cf % 4, cfv)][c.dim] = c.id;

	std::array<int, 3> dim_size = { { charts.nb_charts(0), charts.nb_charts(1), charts.nb_charts(2) } };
	

	// extracting corners relation 
	struct edge {
		const Tetrahedra& _mesh;
		int _c;
		int _cf;
		int _cfv;
		int _chart;
		edge(const Tetrahedra& m, int c, int cf, int cfv, int f) :_mesh(m), _c(c), _cf(cf), _cfv(cfv), _chart(f) {}
		int   _v() {
			return _mesh.facet_vert(_c, _cf, _cfv);
		}
		int   _nv() {
			return _mesh.facet_vert(_c, _cf, (_cfv + 1) % 3);
		}
	};


	std::vector<std::vector<edge>> bnd_chart_vertex_he_graph(m.nverts(), std::vector<edge>());
	{
		std::vector<std::vector<edge>> tmp_vertex_linking(m.nverts(), std::vector<edge>());
		OppositeFacet adj(m);
		FOR(c, m.ncells()) FOR(lcf, 4) if (adj.adjacent[4 * c + lcf] == -1) FOR(cfv, 3)
			tmp_vertex_linking[m.facet_vert(c, lcf, cfv)].push_back(edge(m, c, lcf, cfv, charts.cf2chart[4*c+lcf]));
		FOR(v, m.nverts()) {
			FOR(i, tmp_vertex_linking[v].size()) {
				edge& e1 = tmp_vertex_linking[v][i];
				FOR(j, tmp_vertex_linking[e1._nv()].size()) {
					edge& e2 = tmp_vertex_linking[e1._nv()][j];
					if (e2._nv() == v && e1._chart != e2._chart) {
						bnd_chart_vertex_he_graph[v].push_back(e1);
						break;
					}
				}
			}
		}

	}

	std::vector<bool> is_corners(m.nverts(), false);
	std::vector<std::array<int, 3>> corner_coord(m.nverts(), { -1,-1,-1 });
	FOR(v, m.nverts()) if (v_on_chart[v][0] != -1 && v_on_chart[v][1] != -1 && v_on_chart[v][2] != -1) {
		is_corners[v] = true;
		FOR(dim, 3) corner_coord[v][dim] = charts.dim_rank(v_on_chart[v][dim]);
	}

	std::vector<std::vector<std::vector<std::pair<int, int>>>> pairs_in_X(dim_size[2], std::vector<std::vector<std::pair<int, int>> >(dim_size[1], std::vector<std::pair<int, int>>()));
	FOR(v, m.nverts())  if (is_corners[v]) {
		//FOR(link, bnd_chart_vertex_he_graph[v].size()) {
		for (edge e1 : bnd_chart_vertex_he_graph[v]) {
			std::vector<int> prevent_cycling;
			//edge e1 = bnd_chart_vertex_he_graph[v][link];
			int p_v = e1._v();
			int n_v = e1._nv();
			bool not_looping = true;
			while ((!is_corners[n_v]) && not_looping) {
				int old = p_v;
				p_v = n_v;
				if (bnd_chart_vertex_he_graph[p_v][0]._nv() == old) n_v = bnd_chart_vertex_he_graph[p_v][1]._nv();
				else n_v = bnd_chart_vertex_he_graph[p_v][0]._nv();
				if (prevent_cycling.size() % 30 == 0) for (int i : prevent_cycling) if (i == n_v) {
					Trace::alert("We are looping on a boundary edge. Invalid chart layout. We still try to continue.");
					disp_chart(m, charts[e1._chart], "Invalid_chart");
					not_looping = false;
					break;
				}
				prevent_cycling.push_back(n_v);
			}
			if (corner_coord[v][0] == corner_coord[n_v][0]) continue;


			std::pair<int, int> couple;
			if (corner_coord[v][0] > corner_coord[n_v][0]) {
				couple.first = corner_coord[n_v][0];
				couple.second = corner_coord[v][0];
			}
			else {
				couple.second = corner_coord[n_v][0];
				couple.first = corner_coord[v][0];
			}
			pairs_in_X[corner_coord[v][2]][corner_coord[v][1]].push_back(couple);
		}
	}

	std::vector<std::vector<std::vector<int>>> voxels(dim_size[0] - 1, std::vector<std::vector<int>>(dim_size[1] - 1, std::vector<int>(dim_size[2] - 1, 0)));
	{

		//jordan inspired voxelisation
		std::vector<std::vector<int>> plan_value(dim_size[0] - 1, std::vector<int>(dim_size[1] - 1, 0));
		FOR(Z, dim_size[2] - 1) {

			if (Z > 0) {
				FOR(Y, dim_size[1] - 1) {
					FOR(X, dim_size[0] - 1) {

						voxels[X][Y][Z] = voxels[X][Y][Z - 1];
					}
				}
			}

			std::vector<std::vector<bool>> segments(dim_size[1], std::vector<bool>(dim_size[0], false));
			FOR(Y, dim_size[1]) {
				for (auto couple : pairs_in_X[Z][Y]) {
					for (int X = couple.first; X < couple.second; X++)
						segments[Y][X] = true;
				}

			}

			std::vector<int> line_value(dim_size[0] - 1, 0);
			FOR(Y, dim_size[1] - 1) {
				FOR(X, dim_size[0] - 1) if (segments[Y][X])
					line_value[X]++;
				FOR(X, dim_size[0] - 1) {
					plan_value[X][Y] = line_value[X];
				}
			}
			FOR(Y, dim_size[1] - 1) {
				FOR(X, dim_size[0] - 1) {
					// odd when we are on a chart
					if (plan_value[X][Y] % 2 == 1)
						voxels[X][Y][Z] ++;
				}
			}
		}
	}
	std::vector<std::vector<std::vector<int>>> voxel_vertices_id(dim_size[0], std::vector<std::vector<int>>(dim_size[1], std::vector<int>(dim_size[2], -1)));

	FOR(X, dim_size[0] - 1) {
		FOR(Y, dim_size[1] - 1) {
			FOR(Z, dim_size[2] - 1) {
				if (voxels[X][Y][Z] % 2 == 1) {
					FOR(i, 2) FOR(j, 2) FOR(k, 2)
						voxel_vertices_id[X + i][Y + j][Z + k] = 0;
				}
			}
		}
	}
	std::vector<vec3> voxel_vertices;
	FOR(X, dim_size[0]) {
		FOR(Y, dim_size[1]) {
			FOR(Z, dim_size[2]) {
				if (voxel_vertices_id[X][Y][Z] != -1) {
					voxel_vertices_id[X][Y][Z] = voxel_vertices.size();
					voxel_vertices.push_back(vec3(X, Y, Z));
				}
			}
		}
	}
	blocks.m.points.create_points(voxel_vertices.size());

	FOR(v, blocks.m.nverts()) blocks.block_coord[v] = voxel_vertices[v];
	FOR(v, blocks.m.nverts()) FOR(d, 3) blocks.float_coord[v][d] = charts[charts.ranked(d, (int)voxel_vertices[v][d])].value;
	FOR(v, blocks.m.nverts()) blocks.int_coord[v] = blocks.block_coord[v];

	FOR(v, blocks.m.nverts()) blocks.m.points[v] = blocks.float_coord[v];

	change_points_from_polycuboid_to_mesh(m, polycuboid, *blocks.m.points.data);


	std::vector<std::array<int, 8>> hexes;
	FOR(X, dim_size[0] - 1) {
		FOR(Y, dim_size[1] - 1) {
			FOR(Z, dim_size[2] - 1) {
				if (voxels[X][Y][Z] % 2 == 1) {
					hexes.push_back({
						voxel_vertices_id[X][Y][Z], voxel_vertices_id[X+1][Y][Z], voxel_vertices_id[X][Y+1][Z], voxel_vertices_id[X+1][Y+1][Z],
						voxel_vertices_id[X][Y][Z+1], voxel_vertices_id[X+1][Y][Z+1], voxel_vertices_id[X][Y+1][Z+1], voxel_vertices_id[X+1][Y+1][Z+1]
						});
				}
			}
		}
	}
	blocks.m.create_cells(hexes.size());
	FOR(c, blocks.m.ncells()) FOR(cv, 8) blocks.m.vert(c, cv) = hexes[c][cv];
	FOR(c, blocks.m.ncells()) blocks.cell_volume[c] = blocks.m.util.cell_volume(c);
	
	OppositeFacet hexadj(blocks.m);
	FOR(c, blocks.m.ncells()) FOR(cf, 6) {
		vec3 A = blocks.block_coord[blocks.m.facet_vert(c, cf, 0)];
		vec3 B = blocks.block_coord[blocks.m.facet_vert(c, cf, 1)];
		vec3 C = blocks.block_coord[blocks.m.facet_vert(c, cf, 2)];
		vec3 n = cross(B - A, C - A);
		int d = (std::abs(n * vec3(1, 0, 0)) < 0.1 ? (std::abs(n * vec3(0, 1, 0)) < 0.1 ? 2 : 1) : 0);
		int coord = (int) blocks.block_coord[blocks.m.facet_vert(c, cf, 0)][d];
		blocks.issuing_chart_id[6 * c + cf] = charts.ranked(d, coord);
		if (hexadj.adjacent[6 * c + cf] == -1) blocks.bnd_chart_id[6 * c + cf] = blocks.issuing_chart_id[6 * c + cf];
		else blocks.bnd_chart_id[6 * c + cf] = -1;
	}

	if (!check_mesh_is_connexe(blocks.m)) { Trace::alert("The block decomposition is not connexe. It means that the flagging is bad, we still continue."); }
	if (!check_bnd_mesh_is_variety(blocks.m)) { Trace::alert("The block decomposition boundary is not variety. It means that the flagging is bad, we still continue."); }
}

