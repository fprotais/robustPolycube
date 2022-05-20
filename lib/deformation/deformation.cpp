#include "deformation.h"
#include <utils/trace.h>
#include <OpenNL_psm/OpenNL_psm.h>

#include <utils/elliptic_smoothing.h>

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

void cube_cover(const UM::Tetrahedra& m, const UM::CellFacetAttribute<int>& flag, UM::PointAttribute<UM::vec3>& U, bool constrained) {
	std::cerr << "Initialising polycuboid with Laplacian energy..." << std::endl;
	auto context = nlNewContext();
	DisjointSet ds(m.nverts() * 3);
	FOR(c, m.ncells()) FOR(cf, 4) if (flag[4 * c + cf] != -1) {
		int d = flag[4 * c + cf] / 2;
		FOR(cfv, 3) ds.merge(d * m.nverts() + m.facet_vert(c, cf, cfv), d * m.nverts() + m.facet_vert(c, cf, (cfv + 1) % 3));
	}
	std::vector<int> idmap;
	int nb_var = ds.get_sets_id(idmap);
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_VARIABLES, NLint(nb_var));
	
	nlBegin(NL_SYSTEM);
	if (constrained) {
		FOR(c, m.ncells()) FOR(cf, 4) if (flag[4 * c + cf] != -1) {
			int d = flag[4 * c + cf] / 2;
			FOR(cfv, 3) {
				nlSetVariable(idmap[d * m.nverts() + m.facet_vert(c, cf, cfv)], U[m.facet_vert(c, cf, cfv)][d]);
				nlLockVariable(idmap[d * m.nverts() + m.facet_vert(c, cf, cfv)]);
			}
		}
	}

	nlEnable(NL_VERBOSE);
	nlBegin(NL_MATRIX);
	FOR(c, m.ncells()) {
		int v[4] = { m.vert(c,0) , m.vert(c,1), m.vert(c,2), m.vert(c,3) };
		UM::mat3x3 M = { m.points[v[1]] - m.points[v[0]], m.points[v[2]] - m.points[v[0]], m.points[v[3]] - m.points[v[0]] };
		UM::mat3x3 invM = M.invert();
		double detM = M.det();
		invM = invM.transpose();
		mat<4, 3> grad_coef = { -invM[0] - invM[1] - invM[2], invM[0], invM[1], invM[2] };
		FOR(dim, 3) {
			FOR(dim2, 3) {
				vec3  e(0, 0, 0); e[dim2] = 1;
				nlBegin(NL_ROW);
				FOR(dim_e, 3) FOR(point, 4) {
					nlCoefficient(idmap[dim * m.nverts() + v[point]], e[dim_e] * grad_coef[point][dim_e]);
				}
				if (dim == dim2) nlRightHandSide(1);
				nlEnd(NL_ROW);
			}
		}
	}


	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);
	nlSolve();

	FOR(v, m.nverts()) FOR(dim, 3) U[v][dim] = nlGetVariable(idmap[m.nverts() * dim + v]);

	nlDeleteContext(context);
	std::cerr << " Done.\n";
}



class Polycube_variables : public Elliptic_smoother_Variables {
public:
	Polycube_variables(const Tetrahedra& im, PointAttribute<vec3>& iU, std::vector<int>& iids, std::vector<bool>& ilock_id, int inb_set)
		: m(im)
		, U(iU)
		, ids(iids)
		, lock_id(ilock_id)
		, nb_set(inb_set)
	{
	}
	inline int nb_of_reducted_variables() { return nb_set; }
	inline void set_reducted_values(const std::vector<double>& X) {
		FOR(v, m.nverts()) FOR(d, 3){
			if (!lock_id[ids[3 * v + d]]) U[v][d] = X[ids[3*v+d]];
		}
	}
	inline void get_reducted_values(std::vector<double>& X) {
		FOR(v, m.nverts()) FOR(d, 3) {
			X[ids[3 * v + d]] = U[v][d];
		}
	}
	inline double value(const int i) {
		int cv = i / 3;
		int d = i % 3;
		return U[m.cells[cv]][d];
	}
	inline void add2grad(const int i, const double value, std::vector<double>& G) {
		int cv = i / 3;
		int d = i % 3;
		if (!lock_id[ids[3 * m.cells[cv] + d]]) G[ids[3 * m.cells[cv] + d]] += value;
	}
	inline void add2Hessian(const int i, const int j, const double value, std::vector<sparse_term>& H) {
	} 
private:
	const Tetrahedra& m;
	PointAttribute<vec3>& U;
	std::vector<int>& ids;
	std::vector<bool>& lock_id;
	int nb_set = -1;
};



void correct_param(const UM::Tetrahedra& ref, const UM::CellFacetAttribute<int>& flag, UM::PointAttribute<UM::vec3>& refU, bool constrained) {
	Tetrahedra m;
	*m.points.data = *ref.points.data;
	m.cells = ref.cells;
	UM::PointAttribute<UM::vec3> U(m);
	FOR(v, ref.nverts()) U[v] = refU[v];
	{ // smooth polycube
		DisjointSet ds(3 * m.nverts());
		FOR(c, ref.ncells()) FOR(cf, 4) if (flag[4 * c + cf] != -1) {
			int d = flag[4 * c + cf] / 2;
			FOR(cfv, 3) ds.merge(3 * m.facet_vert(c, cf, cfv) + d, 3 * m.facet_vert(c, cf, (cfv + 1) % 3) + d);
		}
		std::vector<int> ids;
		int nb_sets = ds.get_sets_id(ids);
		std::vector<bool> locks_ids(nb_sets, false);
		if (constrained) {
			FOR(c, ref.ncells()) FOR(cf, 4) if (flag[4 * c + cf] != -1) {
				int d = flag[4 * c + cf] / 2;
				locks_ids[ids[3 * m.facet_vert(c, cf, 0) + d]] = true;
			}
		}

		Polycube_variables var(m, U, ids, locks_ids, nb_sets);

		smoother_options options = _3D_default;
		options.maxiter = 5;
		options.bfgs_maxiter = 1000;
		options.eps_from_theorem = true;
		options.theta = 1. / 4;
		options.bfgs_threshold = 1e-10;
		Elliptic_smoother_3D opt(var, m.ncells(), options);
		opt.start_eps = 1e-2;
		FOR(c, m.ncells()) {
			std::array < UM::vec3, 4> tet_ref;
			FOR(cv, 4) tet_ref[cv] = m.points[m.vert(c, cv)];
			opt.set_tet_ref(c, tet_ref);
		}
		opt.go();

		FOR(v, ref.nverts()) refU[v] = U[v];
	}

}
