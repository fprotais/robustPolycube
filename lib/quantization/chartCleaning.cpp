#include "chartCleaning.h"
#include <algorithm>
#include <utils/trace.h>
#include "tools.h"

#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;
using namespace rb_data_structure;

const std::array<UM::vec3, 6> flag2normal = { UM::vec3(1,0,0), UM::vec3(-1,0,0), UM::vec3(0,1,0), UM::vec3(0,-1,0), UM::vec3(0,0,1), UM::vec3(0,0,-1) };

bool polycuboid_is_valid(const UM::Tetrahedra & m, const UM::Tetrahedra & polycuboid, UM::CellFacetAttribute<int>&cfflag) {
	bool one_is_bad = false;
	CellFacetAttribute<bool> isbad(m, false);
	FOR(cf, m.ncells() * 4) if (cfflag[cf] != -1) {
		int dim = cfflag[cf] / 2;
		int c = cf / 4;
		int lf = cf % 4;
		double x = polycuboid.points[m.facet_vert(c, lf, 0)][dim];
		FOR(i, 2) if (std::abs(x - polycuboid.points[m.facet_vert(c, lf, i+1)][dim]) > 1e-8) {
			isbad[cf] = true;
			one_is_bad = true;
			break;
		}
	}
	if (one_is_bad) {
		Trace::drop_cellfacet_scalar(m, isbad, "NotPlannarFacetInFlagging", false, true);
		return false;
	}
	return true;
}


static void compute_charts(const UM::Triangles & m, const UM::FacetAttribute<int>&flag, const PointAttribute<vec3>& U, Sorted_charts & charts) {

	DisjointSet ds(3 * m.nverts());
	FOR(f, m.nfacets()) FOR(fv, 3) {
		int dim = (flag[f] / 2);
		ds.merge(dim * m.nverts() + m.vert(f, fv), dim * m.nverts() + m.vert(f, (fv + 1) % 3));
	}

	std::vector<int> ids;
	int nb_ids = ds.get_sets_id(ids);
	std::vector<int> bnd_id(nb_ids, -1);


	FOR(f, m.nfacets())  {
		int dim = (flag[f] / 2);
		int rep_v = m.vert(f, 0);
		int id = ids[dim * m.nverts() + rep_v];
		if (bnd_id[id] == -1) {
			bnd_id[id] = charts.size();
			charts.charts.push_back(Chart());
			charts.charts.back().dim = dim;
			charts.charts.back().value = U[rep_v][dim];
			charts.charts.back().tetFacets.push_back(f);
		}
		else {
			charts[bnd_id[id]].tetFacets.push_back(f);
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

	charts.cf2chart.resize(m.nfacets(), -1);
	FOR(c, charts.size()) for (int cf : charts[c].tetFacets) charts.cf2chart[cf] = c;

	FOR(c, charts.size()) charts[c].id = c;

}



inline void mark_corners(const Triangles& m, const Sorted_charts& charts, PointAttribute<bool>& iscorner) {
	std::vector<std::array<bool, 3>> is_on_chart(m.nverts(), { false, false, false });
	FOR(c, charts.size()) for (int f : charts[c].tetFacets) FOR(fv, 3) is_on_chart[m.vert(f, fv)][charts[c].dim] = true;
	FOR(v, m.nverts()) iscorner[v] = is_on_chart[v][0] && is_on_chart[v][1] && is_on_chart[v][2];
}

inline bool make_start_with_a_corner(const Triangles& m, const PointAttribute<bool>& iscorner, const SurfaceConnectivity& fec, std::vector<int>& circular_he) {
	int first_corner = -1;
	FOR(i, circular_he.size()) if (iscorner[fec.from(circular_he[i])]) {
		first_corner = i;
		break;
	}
	if (first_corner == -1) return false;
	std::vector<int> copy = circular_he;
	FOR(i, circular_he.size()) circular_he[i] = copy[(i + first_corner) % circular_he.size()];
	return true;
}

std::vector<int> connexe_chart_componant(const Triangles& m, const SurfaceConnectivity& fec, const Sorted_charts& charts, int source_he) {
	int chart = charts.cf2chart[fec.facet(source_he)];
	DisjointSet ds(m.nfacets());
	FOR(he, m.ncorners()) {
		if (charts.cf2chart[fec.facet(he)] != chart) continue;
		int opp = fec.opposite(he);
		if (charts.cf2chart[fec.facet(opp)] == chart) ds.merge(fec.facet(opp), fec.facet(he));
	}
	std::vector<int> ids;
	ds.get_sets_id(ids);

	std::vector<int> facets;
	int key_id = ids[fec.facet(source_he)];
	FOR(f, m.nfacets()) if (ids[f] == key_id) facets.push_back(f);
	return facets;
}

bool cleanflags(const UM::Tetrahedra& volume, const UM::Tetrahedra& polycuboid, UM::CellFacetAttribute<int>& cfflag) {
	bool flagging_is_locally_valid = false;
	Triangles m;
	std::vector<int> surf2cf;
	OppositeFacet vec(volume);

	m.points.data->assign(volume.points.begin(), volume.points.end());
	FOR(cf, 4 * volume.ncells()) if (vec.adjacent[cf] == -1) {
		int f = m.create_facets(1);
		surf2cf.push_back(cf);
		FOR(fv, 3) m.vert(f, fv) = volume.facet_vert(cf / 4, cf % 4, fv);
	}
	PointAttribute<vec3> U(m);
	FOR(v, m.nverts()) U[v] = polycuboid.points[v];

	m.delete_isolated_vertices();
	FacetAttribute<int> flag(m);
	FOR(f, m.nfacets()) flag[f] = cfflag[surf2cf[f]];

	int MAX_ITER_CORRECTION = 500;
	int iter = 0;

	bool DEFORMATION_LOOKS_BROKEN = false;

	double AVG_EDGE_SIZE = 0;
	FOR(f, m.nfacets()) FOR(fc, 3) {
		AVG_EDGE_SIZE += (U[m.vert(f, fc)] - U[m.vert(f, (fc + 1) % 3)]).norm();
	}
	AVG_EDGE_SIZE /= m.nfacets() * 3;
	

	while (!flagging_is_locally_valid && iter++ < MAX_ITER_CORRECTION) {
		Sorted_charts charts;
		compute_charts(m, flag, U, charts);
		flagging_is_locally_valid = true;

		PointAttribute<bool> iscorner(m);
		mark_corners(m, charts, iscorner);

		std::vector<bool> mark(m.nfacets() * 3, true);
		SurfaceConnectivity fec(m);

		FOR(he, m.ncorners()) {
			int opp = fec.opposite(he);
			if (charts.cf2chart[fec.facet(he)] != charts.cf2chart[fec.facet(opp)]) continue;
			mark[he] = false;
			mark[opp] = false;
		}

		FOR(starting_he, m.ncorners()) if (mark[starting_he]) {
			int actual_chart = charts.cf2chart[fec.facet(starting_he)];
			
			std::vector<int> circular_he;
			int actual_he = starting_he;
			do {
				circular_he.push_back(actual_he);

				int candidate_he = fec.next(actual_he);
				int opp_can = fec.opposite(candidate_he);
				while (charts.cf2chart[fec.facet(opp_can)] == actual_chart) {
					candidate_he = fec.next(opp_can);
					opp_can = fec.opposite(candidate_he);
				}
				actual_he = candidate_he;

			} while (actual_he != starting_he);

			if (!make_start_with_a_corner(m, iscorner, fec, circular_he)) {
				Trace::alert("A chart is isolated into a bigger one, we are merging them.");
				int opp_flag = flag[fec.facet(fec.opposite(circular_he[0]))];
				std::vector<int> badfacets = connexe_chart_componant(m, fec, charts, circular_he[0]);
				for (int f : badfacets) flag[f] = opp_flag;
				flagging_is_locally_valid = false;

				//FacetAttribute<int> here(m, -1);
				//for (int f : badfacets) here[f] = 0;
				//Trace::drop_facet_scalar(m, here, "modified", -1);
				//Trace::drop_facet_scalar(m, flag, "modifiedFlagging");
				break;
 			}

			struct Polycube_edge {
				int c1;
				int c2;
				int dir;
				int opp_chart;
				int opp_dir;
				int opp_dim;
			};
			std::vector<Polycube_edge> pedges;
			
			auto _3rd_chart = [](int chart1, int chart2) -> int {
				int d1 = chart1 / 2;
				int d2 = chart2 / 2;
				int s1 = chart1 % 2;
				int s2 = chart2 % 2;
				if (((d1 + 1) % 3) == d2) {
					int d3 = (d1 + 2) % 3;
					int s3 = (s1 + s2) % 2;
					return 2 * d3 + s3;
				}
				if (((d1 + 2) % 3) == d2) {
					int d3 = (d1 + 1) % 3;
					int s3 = (s1 + s2 + 1) % 2;
					return 2 * d3 + s3;
				}
				return -1;
			};

			for (int he : circular_he) if (iscorner[fec.from(he)]) {
				if (!pedges.empty()){
					pedges.back().c2 = fec.from(he);
				}
				Polycube_edge pedge;
				pedge.c1 = fec.from(he);
				pedge.c2 = -1;
				int opp = fec.opposite(he);
				pedge.opp_chart = charts.cf2chart[fec.facet(opp)];
				pedge.dir = _3rd_chart(flag[fec.facet(he)], flag[fec.facet(opp)]);
				pedge.opp_dim = flag[fec.facet(opp)] / 2;
				pedge.opp_dir = flag[fec.facet(opp)] % 2;
				pedges.push_back(pedge);
			}
			pedges.back().c2 = fec.from(circular_he[0]);
			

			std::array<double, 3> gains = { 0,0,0 };
			std::array<int, 3> neigh_dir = { -2,-2,-2 };

			for (Polycube_edge pedge : pedges) {
				if (neigh_dir[pedge.opp_dim] == -2) neigh_dir[pedge.opp_dim] = pedge.opp_dir;
				else if (neigh_dir[pedge.opp_dim] != pedge.opp_dir) neigh_dir[pedge.opp_dim] = -1;
				
				FOR(d, 3) gains[d] += std::abs(U[pedge.c2][d] - U[pedge.c1][d]);
			}
			int curr_dim = charts[actual_chart].dim;

			FOR(d, 3) if (d != curr_dim) if (gains[d] <= AVG_EDGE_SIZE * 1e-6) {
				int new_dir = neigh_dir[d];
				if (new_dir == -2) {
					Trace::alert("What? trying to continue.");
				}
				else if (new_dir == -1) {
					Trace::alert("The deformation and flagging looks bad. Fixing will introduce invalid flaggings.");
					std::vector<int> badfacets = connexe_chart_componant(m, fec, charts, circular_he[0]);

					FacetAttribute<int> here(m, -1);
					for (int f : badfacets) here[f] = 0;
					Trace::drop_facet_scalar(m, here, "problematic_chart_", -1);
					flagging_is_locally_valid = false;
					DEFORMATION_LOOKS_BROKEN = true;
					break;
				}
				else {
					Trace::alert("A chart is flat in a wrong dimension. Proceeding to some merging to fix that.");
					std::vector<int> badfacets = connexe_chart_componant(m, fec, charts, circular_he[0]);
					for (int f : badfacets) flag[f] = 2 * d + new_dir;
					flagging_is_locally_valid = false;
					break;
				}
			}
			if (!flagging_is_locally_valid) break;

		}

		if (DEFORMATION_LOOKS_BROKEN) break;
	}
	FOR(f, m.nfacets()) cfflag[surf2cf[f]] = flag[f];

	if (!flagging_is_locally_valid) {
		return false;
	}
	else {
		std::cerr << "The flagging looks locally valid." << std::endl;
	}
	return true;
}

