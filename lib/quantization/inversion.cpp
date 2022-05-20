#include "inversion.h"
#include <algorithm>
#include <utils/trace.h>
#include <utils/standard_hexsmoother.h>
#include "tools.h"

#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;
using namespace rb_data_structure;

static void applying_Flm1_with_algorithm_to_extract_hexmesh(const Sorted_charts & charts, const rb_data_structure::Block_decomposition & blocks, Final_mesh& hexmesh) {
	const Hexahedra& m = blocks.m;
	std::vector<std::array<int, 3>> vert_charts(m.nverts(), { -1,-1,-1 });
	FOR(c, m.ncells()) FOR(cf, 6) if (blocks.bnd_chart_id[6 * c + cf] != -1) FOR(cfv, 4)
		vert_charts[m.facet_vert(c, cf, cfv)][charts[blocks.bnd_chart_id[6 * c + cf]].dim] = blocks.bnd_chart_id[6 * c + cf];
	
	// algorithm is here, hidden with a union-find algorithm
	DisjointSet ds(m.nverts());
	FOR(c, m.ncells()) FOR(cf, 6) FOR(cfc, 4) {
		int v0 = m.facet_vert(c, cf, cfc);
		int v1 = m.facet_vert(c, cf, (cfc + 1) % 4);
		double edge_size = (blocks.int_coord[v0] - blocks.int_coord[v1]).norm();
		if (std::round(edge_size) == 0) {
			ds.merge(v0, v1);
		}
	}
	std::vector<int> newid;
	int nb_new_verts = ds.get_sets_id(newid);

	std::vector<std::array<int, 3>> new_vert_charts(nb_new_verts, { -1,-1,-1 });
	FOR(v, m.nverts()) FOR(d, 3) if (vert_charts[v][d] != -1) {
		if ((new_vert_charts[newid[v]][d] != vert_charts[v][d]) && (new_vert_charts[newid[v]][d] != -1)) {
			Trace::alert("Algorithm couldn't go through, invalid quantization. Shouldn't happen under regular use, please investigate.");
			abort();
		}
		new_vert_charts[newid[v]][d] = vert_charts[v][d];
	}
	// end here

	hexmesh.m.points.create_points(nb_new_verts);
	FOR(v, m.nverts()) hexmesh.m.points[newid[v]] = m.points[v];
	FOR(v, m.nverts()) hexmesh.polycube_coord[newid[v]] = blocks.int_coord[v];

	hexmesh.m.create_cells(m.ncells());
	FOR(c, m.ncells()) hexmesh.original_bloc[c] = c;
	FOR(c, m.ncells()) FOR(cv, 8) hexmesh.m.vert(c, cv) = newid[m.vert(c, cv)];

	std::vector<bool> todel(m.ncells(), false);
	FOR(c, hexmesh.m.ncells()) FOR(cf, 6) FOR(cfc, 4) {
		int v0 = hexmesh.m.facet_vert(c, cf, cfc);
		int v1 = hexmesh.m.facet_vert(c, cf, (cfc + 1) % 4);
		// with DisjointSet, it is equivalent to edgesize = 0
		if (v1 == v0) todel[c] = true;
	}
	hexmesh.m.delete_cells(todel);
	hexmesh.m.delete_isolated_vertices();

	FOR(cf, hexmesh.m.ncells() * 6) hexmesh.bnd_chart[cf] = -1;
	FOR(c, hexmesh.m.ncells()) FOR(cf, 6) FOR(d, 3) {
		int chart_id = new_vert_charts[hexmesh.m.facet_vert(c, cf, 0)][d];
		FOR(cfv, 3) {
			int v = hexmesh.m.facet_vert(c, cf, cfv + 1);
			if (new_vert_charts[v][d] != chart_id) {
				chart_id = -1;
				break;
			}
		}
		if (chart_id == -1) continue;
		hexmesh.bnd_chart[6 * c + cf] = chart_id;
		break;
	}

}

inline vec3 baricentric_point(const std::array<vec3, 8>& points, std::array<int, 3> size, const int i, const int j, const int k) {
	vec3 point(0, 0, 0);
	FOR(i, 3) size[i] -= 1;
	const double divisor = (size[0] * size[1] * size[2]);
	point += (size[0] - i) * (size[1] - j) * (size[2] - k) * points[0] / divisor;
	point += (i) * (size[1] - j) * (size[2] - k) * points[1] / divisor;
	point += (size[0] - i) * (j) * (size[2] - k) * points[2] / divisor;
	point += (i) * (j) * (size[2] - k) * points[3] / divisor;
	point += (size[0] - i) * (size[1] - j) * (k)*points[4] / divisor;
	point += (i) * (size[1] - j) * (k)*points[5] / divisor;
	point += (size[0] - i) * (j) * (k)*points[6] / divisor;
	point += (i) * (j) * (k)*points[7] / divisor;
	return point;
}


inline std::array<int, 8> oriented_int_nodes(const Hexahedra& m, const PointAttribute<vec3>& Ui, int c) {
	std::array<int, 3> max = { { 0,0,0 } };
	std::array<int, 3> min = { { -1, -1, -1 } };
	std::array<int, 8> nodes;
	FOR(dim, 3) FOR(cv, 8) {
		if (Ui[m.vert(c, cv)][dim] > max[dim])
			max[dim] = int(Ui[m.vert(c, cv)][dim]);
		else if (Ui[m.vert(c, cv)][dim] < min[dim])
			min[dim] = int(Ui[m.vert(c, cv)][dim]);
	}

	FOR(cv, 8) {
		const vec3 node = Ui[m.vert(c, cv)];
		if (int(node.z) == max[2]) {
			if (int(node.y) == max[1]) {
				if (int(node.x) == max[0])
					nodes[7] = m.vert(c, cv);
				else
					nodes[6] = m.vert(c, cv);
			}
			else {
				if (int(node.x) == max[0])
					nodes[5] = m.vert(c, cv);
				else
					nodes[4] = m.vert(c, cv);
			}
		}
		else {
			if (int(node.y) == max[1]) {
				if (int(node.x) == max[0])
					nodes[3] = m.vert(c, cv);
				else
					nodes[2] = m.vert(c, cv);
			}
			else {
				if (int(node.x) == max[0])
					nodes[1] = m.vert(c, cv);
				else
					nodes[0] = m.vert(c, cv);
			}
		}
	}

	return nodes;
}


static void split_to_unit_hexes(const Final_mesh& coarsehexmesh, Final_mesh& hexmesh) {
	Hexahedra m;
	m.points = coarsehexmesh.m.points;
	m.cells = coarsehexmesh.m.cells;
	
	//assuring we got numering right
	FOR(c, m.ncells()) {
		std::array<int, 8> nodes = oriented_int_nodes(m, coarsehexmesh.polycube_coord, c);
		FOR(cv, 8) m.vert(c, cv) = nodes[cv];
	}

	std::vector<std::array<int, 8>> hexes;
	std::vector<int> original_cell;
	std::vector<vec3> new_pts;
	std::vector<vec3> new_pts_int;

	std::vector<std::array<int, 3>> cell_size(m.ncells());
	std::vector<int> cell_start(m.ncells(), 0);
	auto vert_id = [&cell_size, &cell_start](int c, int k, int j, int i)-> int  {
		return cell_start[c] + k * cell_size[c][0] * cell_size[c][1] + j * cell_size[c][0] + i;
	};
	FOR(c, m.ncells()) {
		std::array<vec3, 8> pts; FOR(cv, 8) pts[cv] = m.points[m.vert(c, cv)];
		//FOR(cv, 8) pts[cv] = coarsehexmesh.polycube_coord[m.vert(c, cv)];
		std::array<vec3, 8> pts_i; FOR(cv, 8) pts_i[cv] = coarsehexmesh.polycube_coord[m.vert(c, cv)];

		cell_size[c] = {
			int(pts_i[1][0] - pts_i[0][0]) + 1,
			int(pts_i[2][1] - pts_i[0][1]) + 1,
			int(pts_i[4][2] - pts_i[0][2]) + 1
		};
		if (c > 0) cell_start[c] = cell_start[c - 1] + cell_size[c - 1][2] * cell_size[c - 1][1] * cell_size[c - 1][0];

		FOR(k, cell_size[c][2]) FOR(j, cell_size[c][1]) FOR(i, cell_size[c][0]) {
			new_pts.push_back(baricentric_point(pts, cell_size[c], i, j, k));
			new_pts_int.push_back(baricentric_point(pts_i, cell_size[c], i, j, k));
		}
		FOR(k, cell_size[c][2] - 1) FOR(j, cell_size[c][1] - 1) FOR(i, cell_size[c][0] - 1) {
			std::array<int, 8> hex;
			FOR(kk, 2) FOR(jj, 2) FOR(ii, 2) hex[4 * kk + 2 * jj + ii] = vert_id(c, k + kk, j + jj, i + ii);
			hexes.push_back(hex);
			original_cell.push_back(coarsehexmesh.original_bloc[c]);
		}
	}

	std::vector<std::array<int, 3>> vert_charts(new_pts.size(), { -1,-1,-1 });
	
	OppositeFacet adj(m);
	DisjointSet ds(new_pts.size());
	FOR(c, m.ncells()) FOR(cf, 6){

		int dim = cf / 2;
		int step = cf % 2;
		if (adj[6 * c + cf] == -1) {
			int ch = coarsehexmesh.bnd_chart[6 * c + cf];
			if (dim == 0) {
				int i = (step == 0 ? 0 : cell_size[c][0] - 1);
				FOR(k, cell_size[c][2]) FOR(j, cell_size[c][1])  vert_charts[vert_id(c, k, j, i)][dim] = ch;
			}
			else if (dim == 1) {
				int j = (step == 0 ? 0 : cell_size[c][1] - 1);
				FOR(k, cell_size[c][2]) FOR(i, cell_size[c][0])  vert_charts[vert_id(c, k, j, i)][dim] = ch;
			}
			else if (dim == 2) {
				int k = (step == 0 ? 0 : cell_size[c][2] - 1);
				FOR(j, cell_size[c][1]) FOR(i, cell_size[c][0])  vert_charts[vert_id(c, k, j, i)][dim] = ch;
			}
		}
		else {
			int opp_c = adj[6 * c + cf] / 6;
			if (dim == 0) {
				int i = 0;
				int opp_i = cell_size[opp_c][0] - 1;
				if (step == 1) {
					i = cell_size[c][0] - 1;
					opp_i = 0;
				}
				FOR(k, cell_size[c][2]) FOR(j, cell_size[c][1]) 
					ds.merge(vert_id(c, k, j, i), vert_id(opp_c, k, j, opp_i));
			}
			else if (dim == 1) {
				int j = 0;
				int opp_j = cell_size[opp_c][1] - 1;
				if (step == 1) {
					j = cell_size[c][1] - 1;
					opp_j = 0;
				}
				FOR(k, cell_size[c][2]) FOR(i, cell_size[c][0])
					ds.merge(vert_id(c, k, j, i), vert_id(opp_c, k, opp_j, i));
			}
			else if (dim == 2) {
				int k = 0;
				int opp_k = cell_size[opp_c][2] - 1;
				if (step == 1) {
					k = cell_size[c][2] - 1;
					opp_k = 0;
				}
				FOR(j, cell_size[c][1]) FOR(i, cell_size[c][0])
					ds.merge(vert_id(c, k, j, i), vert_id(opp_c, opp_k, j, i));
			}
		}
	}
	std::vector<int> new_vert_id;
	int nb_new_vert = ds.get_sets_id(new_vert_id);

	hexmesh.m.points.create_points(nb_new_vert);
	FOR(v, new_pts.size()) hexmesh.m.points[new_vert_id[v]] = new_pts[v];
	FOR(v, new_pts.size()) hexmesh.polycube_coord[new_vert_id[v]] = new_pts_int[v];
	hexmesh.m.create_cells(hexes.size());
	FOR(c, hexmesh.m.ncells()) FOR(cv, 8) hexmesh.m.vert(c, cv) = new_vert_id[hexes[c][cv]];
	FOR(c, hexmesh.m.ncells()) hexmesh.original_bloc[c] = original_cell[c];


	std::vector<std::array<int, 3>> swaped_chart(nb_new_vert, {-1,-1,-1});
	FOR(v, new_pts.size()) FOR(d,3) if (swaped_chart[new_vert_id[v]][d] == -1) swaped_chart[new_vert_id[v]][d] = vert_charts[v][d];

	FOR(cf, hexmesh.m.ncells() * 6) hexmesh.bnd_chart[cf] = -1;
	FOR(c, hexmesh.m.ncells()) FOR(cf, 6) FOR(d, 3) {
		int chart_id = swaped_chart[hexmesh.m.facet_vert(c, cf, 0)][d];
		FOR(cfv, 3) {
			int v = hexmesh.m.facet_vert(c, cf, cfv + 1);
			if (swaped_chart[v][d] != chart_id) {
				chart_id = -1;
				break;
			}
		}
		if (chart_id == -1) continue;
		hexmesh.bnd_chart[6 * c + cf] = chart_id;
		break;
	}
}

static void smooth_in_polycuboid(const Sorted_charts& charts, Final_mesh& hexmesh) {
	std::vector<bool> locks(hexmesh.m.nverts()*3, false);
	FOR(c, hexmesh.m.ncells()) FOR(cf, 6) {
		int ch = hexmesh.bnd_chart[6 * c + cf];
		if (ch == -1) continue;
		int dim = charts[ch].dim;
		FOR(cfv, 4) hexmesh.m.points[hexmesh.m.facet_vert(c, cf, cfv)][dim] = charts[ch].value;
		FOR(cfv, 4) locks[3 * hexmesh.m.facet_vert(c, cf, cfv) + dim] = true;
	}
	smooth_hex_mesh(hexmesh.m, locks);
}

void inverse(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, const Sorted_charts& charts, const Block_decomposition& blocks, Final_mesh& hexmesh) {
	Final_mesh coarsehexmesh;
	
	applying_Flm1_with_algorithm_to_extract_hexmesh(charts, blocks, coarsehexmesh);

	Trace::drop_cells_scalar(coarsehexmesh.m, coarsehexmesh.original_bloc, "coarsehexmesh");
	Trace::drop_cellfacet_scalar(coarsehexmesh.m, coarsehexmesh.bnd_chart, "coarsehexmesh_charts", -1, true);
	split_to_unit_hexes(coarsehexmesh, hexmesh);

	smooth_in_polycuboid(charts, hexmesh);
	change_points_from_polycuboid_to_mesh(m, polycuboid, *hexmesh.m.points.data);
}