#include "CADaware_hexsmoothing.h"
#include "tools.h"

#include <utils/projection.h>
#include <utils/trace.h>


#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;

void decompose_bnd(const Triangles & m, std::vector<std::vector<int>>&composants, FacetAttribute<int>&facet_id) {
	SurfaceConnectivity fec(m);
	DisjointSet ds(m.nfacets());

	FOR(he, m.nfacets() * 3) {
		int opp_he = fec.opposite(he);
		if (opp_he == -1) continue;
		int f1 = fec.facet(he);
		int f2 = fec.facet(opp_he);
		vec3 n1 = m.util.normal(f1);
		vec3 n2 = m.util.normal(f2);
		if (n1 * n2 > 0.9)   ds.merge(f1, f2);
	}
	std::vector<int> ids;
	int nb_set = ds.get_sets_id(ids);
	composants = {};
	composants.resize(nb_set);
	FOR(f, m.nfacets()) {
		composants[ids[f]].push_back(f);
		facet_id[f] = ids[f];
	}
}

void get_submesh(const Triangles & main, const std::vector<int>&composant, Triangles & sub) {
	sub = Triangles();
	*sub.points.data = *main.points.data;
	sub.create_facets(composant.size());
	FOR(i, sub.nfacets()) FOR(fv, 3) sub.vert(i, fv) = main.vert(composant[i], fv);
}


void project(const Triangles & bnd, Hexahedra & hexm, PointAttribute<int>&bnd_points) {
	std::vector<BBox3> boxes(bnd.nfacets());
	FOR(f, bnd.nfacets()) FOR(fv, 3) boxes[f].add(bnd.points[bnd.vert(f, fv)]);

	HBoxes looker(boxes);
	FOR(v, hexm.nverts()) if (bnd_points[v]) {
		BBox3 box;
		box.add(hexm.points[v]);
		std::vector<int> possible_tri;
		double dil = 0.1;
		do {
			box.dilate(dil);
			looker.intersect(box, possible_tri);
			double best_d = 1E100;
			vec3 best_points = hexm.points[v];

			for (int t : possible_tri) {
				std::array<vec3, 3> ABC;
				FOR(fv, 3) ABC[fv] = bnd.points[bnd.vert(t, fv)];
				vec3 p; std::array<double, 3> l;
				double here = point_triangle_squared_distance(hexm.points[v], ABC, p, l);
				if (here < best_d) {
					best_d = here;
					best_points = p;
				}
			}
			hexm.points[v] = best_points;
			dil *= 2;
		} while (possible_tri.empty());
	}
}

void get_bnd_componant(const Triangles & bnd, const Hexahedra & hexm, const PointAttribute<int>&bnd_points, const std::vector<std::vector<int>>&composants, std::vector<std::vector<int>>&proj, double threshold = 1E-6) {
	FOR(c, composants.size()) {
		Triangles m;
		get_submesh(bnd, composants[c], m);
		std::vector<BBox3> boxes(m.nfacets());
		FOR(f, m.nfacets()) FOR(fv, 3) boxes[f].add(m.points[m.vert(f, fv)]);
		HBoxes looker(boxes);

		FOR(v, hexm.nverts()) if (bnd_points[v]) {
			BBox3 box;
			box.add(hexm.points[v]);
			std::vector<int> possible_tri;
			double dil = threshold;
			box.dilate(dil);
			looker.intersect(box, possible_tri);

			for (int t : possible_tri) {
				std::array<vec3, 3> ABC;
				FOR(fv, 3) ABC[fv] = m.points[m.vert(t, fv)];
				vec3 p; std::array<double, 3> l;
				double d = point_triangle_squared_distance(hexm.points[v], ABC, p, l);
				if (d < threshold) {
					proj[v].push_back(c);
					break;
				}
			}
		}
	}
}

constexpr double w1 = 1. / 27;
constexpr double w2 = 1. / (2 * 27);
constexpr double w3 = 1. / (4 * 27);
constexpr double w4 = 1. / (8 * 27);

struct quadrature {
	double weight;
	int points[3][2];
};

constexpr quadrature HEX_QUADRATURE[64] = {
	// cat I
	{w1, {{0,1},{0,2},{0,4}}},
	{w1, {{1,3},{1,0},{1,5}}},
	{w1, {{2,0},{2,3},{2,6}}},
	{w1, {{3,2},{3,1},{3,7}}},
	{w1, {{4,6},{4,5},{4,0}}},
	{w1, {{5,4},{5,7},{5,1}}},
	{w1, {{6,7},{6,4},{6,2}}},
	{w1, {{7,5},{7,6},{7,3}}},

	// cat II
	{w2, {{4,6},{0,4},{0,1}}},
	{w2, {{6,7},{2,6},{2,0}}},
	{w2, {{7,5},{3,7},{3,2}}},
	{w2, {{5,4},{1,5},{1,3}}},

	{w2, {{6,2},{4,6},{4,5}}},
	{w2, {{2,3},{0,2},{0,4}}},
	{w2, {{3,7},{1,3},{1,0}}},
	{w2, {{7,6},{5,7},{5,1}}},

	{w2, {{0,4},{2,0},{2,3}}},
	{w2, {{4,5},{6,4},{6,2}}},
	{w2, {{5,1},{7,5},{7,6}}},
	{w2, {{1,0},{3,1},{3,7}}},

	{w2, {{0,2},{1,0},{1,5}}},
	{w2, {{2,6},{3,2},{3,1}}},
	{w2, {{6,4},{7,6},{7,3}}},
	{w2, {{4,0},{5,4},{5,7}}},

	{w2, {{5,7},{4,5},{4,0}}},
	{w2, {{7,3},{6,7},{6,4}}},
	{w2, {{3,1},{2,3},{2,6}}},
	{w2, {{1,5},{0,1},{0,2}}},

	{w2, {{1,3},{5,1},{5,4}}},
	{w2, {{3,2},{7,3},{7,5}}},
	{w2, {{2,0},{6,2},{6,7}}},
	{w2, {{0,1},{4,0},{4,6}}},

	// cat III
	{w3, {{1,0},{0,4},{5,7}}},
	{w3, {{5,1},{1,0},{4,6}}},
	{w3, {{4,5},{5,1},{0,2}}},
	{w3, {{0,4},{4,5},{1,3}}},

	{w3, {{3,1},{1,5},{7,6}}},
	{w3, {{3,1},{3,7},{5,4}}},
	{w3, {{7,3},{7,5},{1,0}}},
	{w3, {{5,7},{5,1},{3,2}}},

	{w3, {{0,2},{0,4},{6,7}}},
	{w3, {{0,2},{2,6},{4,5}}},
	{w3, {{6,4},{6,2},{0,1}}},
	{w3, {{4,0},{4,6},{2,3}}},

	{w3, {{6,2},{2,3},{7,5}}},
	{w3, {{2,3},{3,7},{6,4}}},
	{w3, {{3,7},{7,6},{2,0}}},
	{w3, {{7,6},{6,2},{3,1}}},

	{w3, {{0,1},{0,2},{3,7}}},
	{w3, {{0,1},{1,3},{2,6}}},
	{w3, {{1,3},{3,2},{0,4}}},
	{w3, {{2,0},{2,3},{1,5}}},

	{w3, {{5,4},{4,6},{7,3}}},
	{w3, {{5,4},{5,7},{6,2}}},
	{w3, {{7,5},{7,6},{4,0}}},
	{w3, {{4,6},{6,7},{5,1}}},

	// cat IV
	{w4, {{1,0},{2,6},{5,7}}},
	{w4, {{3,1},{0,4},{7,6}}},
	{w4, {{2,3},{1,5},{6,4}}},
	{w4, {{0,2},{3,7},{4,5}}},

	{w4, {{1,0},{3,7},{4,6}}},
	{w4, {{3,1},{2,6},{5,4}}},
	{w4, {{2,3},{0,4},{7,5}}},
	{w4, {{0,2},{1,5},{6,7}}},
};
vec3 unit_cube[8] = { vec3(0,0,0), vec3(1,0,0),vec3(0,1,0), vec3(1,1,0), vec3(0,0,1), vec3(1,0,1), vec3(0,1,1), vec3(1,1,1) };
inline double chi(double eps, double det) {
	if (det > 0)
		return (det + std::sqrt(eps * eps + det * det)) * .5;
	return .5 * eps * eps / (std::sqrt(eps * eps + det * det) - det);
}

inline double chi_deriv(double eps, double det) {
	return .5 + det / (2. * std::sqrt(eps * eps + det * det));
}
inline UM::mat3x3 dual_basis(const UM::mat3x3 & J) {
	return
	{
		{{
			 J[1].y * J[2].z - J[1].z * J[2].y,
			 J[1].z * J[2].x - J[1].x * J[2].z,
			 J[1].x * J[2].y - J[1].y * J[2].x
		 },
		{
			J[0].z * J[2].y - J[0].y * J[2].z,
			J[0].x * J[2].z - J[0].z * J[2].x,
			J[0].y * J[2].x - J[0].x * J[2].y
		},
		{
			J[0].y * J[1].z - J[0].z * J[1].y,
			J[0].z * J[1].x - J[0].x * J[1].z,
			J[0].x * J[1].y - J[0].y * J[1].x
		}}
	};
}
struct lisseur {
	lisseur(Hexahedra& m, std::vector<bool>& locks) : m_(m), locks_(locks), cell_mindet(m_) {
		update_cell_locking();
		double E = compute_elliptic_energy();
		std::cerr << "Lisseur startup, energy: " << E << " eps: " << glob_eps << " detmin: " << mindet << " ninv: " << nb_inverted << std::endl;
		bnd_proj.resize(m_.nverts());
		glob_eps = std::sqrt(min_glob_eps * min_glob_eps + 0.04 * std::min(mindet, 0.) * std::min(mindet, 0.));

	}
	Hexahedra& m_;
	std::vector<bool>& locks_;
	CellAttribute<double> cell_mindet;
	std::vector<bool> locked_cell_;
	std::vector<bool> quadrature2skip_;

	double bnd_w = 1;

	double glob_eps = 1e-6;
	int maxiter = 5;
	double static_threshold = 1E-8;
	double mindet = 1E100;
	int nb_inverted = 0;

	bool debug = true;

	int NB_QUADRATURE = 8;

	double min_glob_eps = 1e-9;

	double lbfgs_threshold = 1E-16;
	int lbfgs_maxiter = 20;

	void update_cell_locking() {
		//CellAttribute<int> c_is_in(m_, 0);
		quadrature2skip_.resize(m_.ncells() * NB_QUADRATURE);
		locked_cell_.resize(m_.ncells(), true);
		FOR(c, m_.ncells()) FOR(i, NB_QUADRATURE) {
			bool is_in = false;
			FOR(d, 3) if (!locks_[3 * m_.vert(c, HEX_QUADRATURE[i].points[d][1]) + d]) is_in = true;
			FOR(d, 3) if (!locks_[3 * m_.vert(c, HEX_QUADRATURE[i].points[d][0]) + d]) is_in = true;
			quadrature2skip_[c * NB_QUADRATURE + i] = !is_in;
			//if (!quadrature2skip_[c * NB_QUADRATURE + i]) c_is_in[c]++;
			if (!quadrature2skip_[c * NB_QUADRATURE + i]) locked_cell_[c] = false;
		}
		//Trace::drop_cells_scalar(m_, c_is_in, "is_in");
		//exit(0);
	}

	double compute_elliptic_energy() {
		double F = 0;
		nb_inverted = 0;
		mindet = 1E100;

		FOR(c, m_.ncells()) if (!locked_cell_[c]) {
			std::array<vec3, 8> hex;
			FOR(i, 8) hex[i] = m_.points[m_.vert(c, i)];
			cell_mindet[c] = 1E100;
			double loc_eps = std::sqrt(glob_eps * glob_eps + 0.04 * std::min(cell_mindet[c], 0.) * std::min(cell_mindet[c], 0.));


			FOR(q, NB_QUADRATURE) if (!quadrature2skip_[c * NB_QUADRATURE + q]) {
				mat3x3 	J = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
				FOR(d, 3) J[d] = hex[HEX_QUADRATURE[q].points[d][1]] - hex[HEX_QUADRATURE[q].points[d][0]];
				double trJJ = J[0] * J[0] + J[1] * J[1] + J[2] * J[2];
				double detJ = J.det();

				double c1 = chi(loc_eps, detJ);
				double c2 = std::pow(c1, 2. / 3.);

				double f = trJJ / c2;
				F += HEX_QUADRATURE[q].weight * f;

				mindet = std::min(mindet, detJ);
				cell_mindet[c] = std::min(cell_mindet[c], detJ);

			}
			nb_inverted += (cell_mindet[c] <= 0);

		}
		return F;
	}

	double compute_bnd_energy() {
		double E = 0;
		FOR(v, m_.nverts()) for (int c : bnd_proj[v]) {
			//double d = rbfs[c].eval(m_.points[v]);
			vec3 proj;
			double d = get_closet_on_componant(c, m_.points[v], proj);
			E += d;
		}
		return bnd_w * E;
	}

	double compute_energy_with_grad(const std::vector<double>& coord, std::vector<double>& grad) {
		double F = 0;

		FOR(c, m_.ncells()) if (!locked_cell_[c]) {
			std::array<vec3, 8> hex;
			double loc_eps = std::sqrt(glob_eps * glob_eps + 0.04 * std::min(cell_mindet[c], 0.) * std::min(cell_mindet[c], 0.));

			FOR(i, 8) FOR(d, 3) hex[i][d] = coord[3 * m_.vert(c, i) + d];

			FOR(q, NB_QUADRATURE) if (!quadrature2skip_[c * NB_QUADRATURE + q]) {
				mat3x3 	J = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
				FOR(d, 3) J[d] = hex[HEX_QUADRATURE[q].points[d][1]] - hex[HEX_QUADRATURE[q].points[d][0]];
				double trJtJ = J[0] * J[0] + J[1] * J[1] + J[2] * J[2];
				double detJ = J.det();

				double c1 = chi(loc_eps, detJ);
				double c2 = std::pow(c1, 2. / 3.);
				double c3 = chi_deriv(loc_eps, detJ);

				mat3x3 K = dual_basis(J);

				double f = trJtJ / c2;
				F += HEX_QUADRATURE[q].weight * f;
				mat3x3 df_dJ = J * (2. / c2) - K * ((2. * f * c3) / (3. * c1));

				FOR(d, 3) FOR(dimgrad, 3) {
					int v1 = 3 * m_.vert(c, HEX_QUADRATURE[q].points[d][1]) + dimgrad;
					int v2 = 3 * m_.vert(c, HEX_QUADRATURE[q].points[d][0]) + dimgrad;
					if (!locks_[v1]) grad[v1] += HEX_QUADRATURE[q].weight * df_dJ[d][dimgrad];
					if (!locks_[v2]) grad[v2] -= HEX_QUADRATURE[q].weight * df_dJ[d][dimgrad];
				}
			}
		}
		double bnd_E = 0;

		FOR(v, m_.nverts()) for (int c : bnd_proj[v]) {
			vec3 p = { coord[3 * v + 0], coord[3 * v + 1], coord[3 * v + 2] };
			vec3 proj;
			double d = get_closet_on_componant(c, p, proj);
			vec3 g = proj - p;
			bnd_E += d;
			vec3 dEdv = -bnd_w * 2 * g;
			FOR(dim, 3) {
				int id = 3 * v + dim;
				if (!locks_[id]) grad[id] += dEdv[dim];
			}
		}

		return F + bnd_w * bnd_E;
	}
	void callback(int iter) {
		if (debug) Trace::drop_cells_scalar(m_, cell_mindet, "det_iter_" + std::to_string(iter) + "_");
	}

	bool update_coord() {
		std::vector<double> newcoord(3 * m_.nverts());
		FOR(v, m_.nverts()) FOR(d, 3) newcoord[3 * v + d] = m_.points[v][d];

		bool res = run_lbfgs(newcoord);


		FOR(v, m_.nverts()) FOR(d, 3) m_.points[v][d] = newcoord[3 * v + d];
		return res;
	}

	void update_eps(double sigma) {
		double mu = (1 - sigma) * chi(glob_eps, mindet);
		if (mindet < mu)
			glob_eps = 2 * std::sqrt(mu * (mu - mindet));
		else glob_eps = min_glob_eps;
	}

	void go() {
		update_cell_locking();
		double E_prev = compute_elliptic_energy(); // compute mindet
		update_eps(0.5);

		if (debug) std::cerr << "Init energy: " << E_prev << " glob eps: " << glob_eps << " detmin: " << mindet << " ninv: " << nb_inverted << std::endl;
		callback(0);
		FOR(iter, maxiter) {
			if (debug) 	std::cerr << "iteration #" << iter << std::endl;
			E_prev = compute_elliptic_energy();
			double bnde_prev = compute_bnd_energy();
			bool res = update_coord();

			double E = compute_elliptic_energy();
			double bndE = compute_bnd_energy();

			if (debug) std::cerr << " E : " << E_prev << "  --->  " << E << " glob eps: " << glob_eps << " detmin: " << mindet << " ninv: " << nb_inverted << " | |E_prev - E| / E = " << std::abs(E_prev - E) / E << std::endl;
			if (debug) std::cerr << "Bnd : " << bnde_prev << " ---> " << bndE << std::endl;
			double sigma = std::max(1. - E / E_prev, 0.1);
			update_eps(sigma);
			callback(iter + 1);
			if (res && glob_eps == min_glob_eps && mindet > 0) break;
		}
		if (debug) std::cerr << "Final state => E : " << compute_elliptic_energy() << " eps: " << glob_eps << " detmin: " << mindet << " ninv: " << nb_inverted << std::endl;

	}



	bool run_lbfgs(std::vector<double>& coord) {
		const STLBFGS::Optimizer::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
			std::fill(G.begin(), G.end(), 0);
			F = compute_energy_with_grad(X, G);
		};
		STLBFGS::Optimizer opt{ func };
		opt.gtol = lbfgs_threshold;
		opt.ftol = lbfgs_threshold;
		opt.maxiter = lbfgs_maxiter;
		opt.verbose = debug;
		opt.invH.history_depth = 10;
		opt.run(coord);
		return true;
	}

	std::vector<HBoxes> componant_lookers;
	std::vector<std::vector<int>> bnd_proj;
	std::vector<std::vector<int>> componants;
	const Triangles* bnd = nullptr;
	void set_bnd_ref(const Triangles& m, const std::vector<std::vector<int>>& composants, const std::vector<std::vector<int>>& vert_componant) {
		bnd = &m;
		bnd_proj = vert_componant;
		componants = composants;
		FOR(c, composants.size()) {
			std::vector<BBox3> boxes(composants[c].size());
			FOR(f, composants[c].size()) FOR(fv, 3) boxes[f].add(bnd->points[bnd->vert(composants[c][f], fv)]);
			componant_lookers.push_back(HBoxes(boxes));
		}
	}

	double get_closet_on_componant(const int c, const vec3 p, vec3& closet_point) {

		BBox3 box;
		box.add(p);
		std::vector<int> possible_tri;
		double dil = 1e-6;
		double best_d = 1E100;
		do {
			box.dilate(dil);
			componant_lookers[c].intersect(box, possible_tri);
			closet_point = p;

			for (int t : possible_tri) {
				std::array<vec3, 3> ABC;
				FOR(fv, 3) ABC[fv] = bnd->points[bnd->vert(componants[c][t], fv)];
				vec3 newp; std::array<double, 3> l;
				double here = point_triangle_squared_distance(p, ABC, newp, l);
				if (here < best_d) {
					best_d = here;
					closet_point = newp;
				}
			}

			dil *= 2;
			if (dil > 1E100) break;
		} while (possible_tri.empty());
		return best_d;
	}

};

void smooth(UM::Hexahedra & hexm, const UM::Triangles & bnd, int debug) {
	const Triangles& m = bnd;
	std::vector<std::vector<int>> composants;
	FacetAttribute<int> facet_id(m);
	decompose_bnd(m, composants, facet_id);
	if (debug > 1)Trace::drop_facet_scalar(m, facet_id, "ids");

	FOR(c, composants.size()) {
		Triangles sub;
		get_submesh(m, composants[c], sub);
		if (debug > 1) Trace::drop_surface(sub, "comp", {});
	}


	PointAttribute<int> bndpoints(hexm, false);
	OppositeFacet vec(hexm);
	FOR(c, hexm.ncells()) FOR(cf, 6) if (vec.adjacent[6 * c + cf] == -1) FOR(cfv, 4) bndpoints[hexm.facet_vert(c, cf, cfv)] = true;


	project(m, hexm, bndpoints);
	std::vector<std::vector<int>> bnd_proj(hexm.nverts());


	get_bnd_componant(m, hexm, bndpoints, composants, bnd_proj, 5E-4 * avgedgesize(hexm));
	PointAttribute<int> proj_disp(hexm);

	FOR(v, hexm.nverts()) proj_disp[v] = bnd_proj[v].size();

	if (debug) Trace::drop_volume_points_scalar(hexm, proj_disp, "number_of_componant", -2);

	FOR(v, hexm.nverts()) if (bnd_proj[v].empty()) proj_disp[v] = -1;
	else proj_disp[v] = bnd_proj[v][0];

	if (debug > 1) Trace::drop_volume_points_scalar(hexm, proj_disp, "bnd_componant", -2);

	std::vector<bool> locks(hexm.nverts() * 3, false);
	{
		PointAttribute<int> distfrombnd(hexm, 10);
		FOR(v, hexm.nverts()) if (bndpoints[v]) distfrombnd[v] = 0;
		FOR(i, 10) FOR(h, hexm.ncells()) FOR(hf, 6) FOR(hfc, 4) FOR(j, 2) {
			int shift = (j + 1) % 2;
			int v1 = hexm.facet_vert(h, hf, (hfc + j) % 4);
			int v2 = hexm.facet_vert(h, hf, (hfc + shift) % 4);
			if (distfrombnd[v1] == i && distfrombnd[v2] > i) distfrombnd[v2] = i + 1;
		}
		FOR(v, hexm.nverts()) if (distfrombnd[v] > 2) FOR(d, 3) locks[3 * v + d] = true;
		if (debug > 1) Trace::drop_volume_points_scalar(hexm, distfrombnd, "dist");
	}
	//FOR(v, hexm.nverts()) if (!bndpoints[v]) FOR(d, 3) locks[3 * v + d] = true;

	lisseur opt(hexm, locks);
	opt.set_bnd_ref(m, composants, bnd_proj);
	opt.lbfgs_maxiter = 30;
	opt.bnd_w = 10000.;
	opt.maxiter = 15;
	opt.min_glob_eps = 1e-14;
	opt.debug = debug > 1;
	std::cerr << "Smoothing init...";
	opt.go();
	if (debug) Trace::drop_cells_scalar(opt.m_, opt.cell_mindet, "det_iter_" + std::to_string(0) + "_");
	std::cerr << "ok. Min det: " << opt.mindet << std::endl;
	FOR(i, 5) {
		opt.glob_eps = 1e-10;
		//opt.bnd_w *= 2;
		opt.lbfgs_maxiter = 20;
		opt.maxiter = 6;
		opt.debug = debug > 1;
		std::cerr << "Smoothing iter " << i << "...";
		opt.go();
		std::cerr << "ok. Min det: " << opt.mindet << std::endl;
		if (debug) Trace::drop_cells_scalar(opt.m_, opt.cell_mindet, "det_iter_" + std::to_string(i + 1) + "_");

	}

}