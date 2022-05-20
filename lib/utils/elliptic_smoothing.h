#pragma once

#include <ultimaille/algebra/mat.h>

#include <vector>
#include <array>
#include <assert.h>

class  Elliptic_smoother_Variables {
public:
	virtual int nb_of_reducted_variables() = 0;
	virtual void set_reducted_values(const std::vector<double>& X) = 0;
	virtual void get_reducted_values(std::vector<double>& X) = 0;

	virtual double value(const int i) = 0;
	virtual void add2grad(const int i, const double value, std::vector<double>& G) = 0;
	// Do not need to be unique in H
	typedef std::pair<std::array<int, 2>, double> sparse_term;
	virtual void add2Hessian(const int i, const int j, const double value, std::vector<sparse_term>& H) = 0;
};



struct smoother_options {
	smoother_options(double _theta, int _maxiter, double _bfgs_threshold, int _bfgs_maxiter, int _debug, bool eps_from_theorem_, bool stopping_when_static_, double static_threshold_, bool _barrier, bool _use_newton)
		: theta(_theta)
		, maxiter(_maxiter)
		, bfgs_threshold(_bfgs_threshold)
		, bfgs_maxiter(_bfgs_maxiter)
		, debug(_debug)
		, eps_from_theorem(eps_from_theorem_)
		, stopping_when_static(stopping_when_static_)
		, static_threshold(static_threshold_)
		, barrier(_barrier)
		, use_newton(_use_newton)
	{}
	double theta = 1./128;
	int maxiter = 10000;
	double bfgs_threshold = .1;
	int bfgs_maxiter = 30000;
	int debug = 1;
	bool eps_from_theorem = false;
	bool stopping_when_static = false;
	double static_threshold = 1e-5;
	bool barrier = false;
	bool use_newton = false;

};

const smoother_options _2D_default(1. / 128, 500, 1e-9, 30000, 1, true, false, 1e-9, false, false);
const smoother_options _3D_default(1. / 128., 500, 1e-9, 300, 1, true, false, 1e-9, false, false);


class Elliptic_smoother_2D {
public:
	Elliptic_smoother_2D(Elliptic_smoother_Variables& var, const int nb_tri, const smoother_options& options = _2D_default)
		: var_(var)
		, N_tri_(nb_tri)
		, options_(options)
	{
		init();
	}
	// v1x, v1y, v2x, v2y, v3x, v3y
	void set_triangle_ref(const int t, const std::array<double, 6> &ref);
	void set_triangle_ref(const std::vector<std::array<double, 6>> &refs);
	void set_area(const int t, const double area);
	bool go();


public:
	virtual inline void energy(const std::vector<double>& X, double& F, std::vector<double>& G) {
		standard_elliptic_energy(X, F, G);
	}
	
	virtual inline void compute_energy_hessian(const std::vector<double>& X, double& F, std::vector<double>& G, std::vector<Elliptic_smoother_Variables::sparse_term>& H) {
		standard_elliptic_energy_w_hessian(X, F, G, H);
	}

	void standard_elliptic_energy(const std::vector<double>& X, double& F, std::vector<double>& G);
	void standard_elliptic_energy_w_hessian(const std::vector<double>& X, double& F, std::vector<double>& G, std::vector<Elliptic_smoother_Variables::sparse_term>& H);
	void init();


	typedef std::array<double, 6> tri_grad;
	typedef std::array<double, 4> mat22;
	virtual double evaluate_energy();
	void evaluate_jacobian();

	bool run_lbfgs(std::vector<double>& X);
	bool run_newton(std::vector<double>& X);

	// inputs
	Elliptic_smoother_Variables& var_;
	int N_tri_;
	smoother_options options_;
	std::vector<tri_grad> refs_grad_;


	// intern variables 
	std::vector<mat22> J_; // per-tet Jacobian matrix = [[JX.x JX.y, JX.z], [JY.x, JY.y, JY.z], [JZ.x, JZ.y, JZ.z]]
	std::vector<mat22> K_; // per-tet dual basis: det J = dot J[i] * K[i]
	std::vector<double> det_; // per-tet determinant of the Jacobian matrix
	std::vector<double> area_; // per-tet determinant of the Jacobian matrix
	double eps_;       // regularization parameter, depends on min(jacobian)

	double detmin_;    // min(jacobian) over all tetrahedra
	int ninverted_; // number of inverted tetrahedra

	double start_eps = 1e-3;

};


class Elliptic_smoother_3D {
public:
	Elliptic_smoother_3D(Elliptic_smoother_Variables& var, const int nb_tets, const smoother_options& options = _3D_default)
		: var_(var)
		, N_tets_(nb_tets)
		, options_(options)
	{
		init();
	}
	void set_tet_ref(const int t, const std::array < UM::vec3, 4> tet_ref);
	void set_tet_ref(const std::vector<std::array < UM::vec3, 4> > tets_ref);
	bool go();


public:
	virtual inline void energy(const std::vector<double>& X, double& F, std::vector<double>& G) {
		standard_elliptic_energy(X, F, G);
	}
	virtual inline void iter_call_back(int iter_nb) {
		(void) iter_nb;
	}
	virtual inline void compute_energy_hessian(const std::vector<double>& X, double& F, std::vector<double>& G, std::vector<Elliptic_smoother_Variables::sparse_term>& H) {
		standard_elliptic_energy_w_hessian(X, F, G, H);
	}

	void standard_elliptic_energy(const std::vector<double>& X, double& F, std::vector<double>& G);
	void standard_elliptic_energy_w_hessian(const std::vector<double>& X, double& F, std::vector<double>& G, std::vector<Elliptic_smoother_Variables::sparse_term>& H);
	void init();


	double evaluate_energy();
	void evaluate_jacobian();

	bool run_lbfgs(std::vector<double>& X);
	bool run_newton(std::vector<double>& X);

	// inputs
	Elliptic_smoother_Variables& var_;
	int N_tets_;
	const smoother_options options_;
	std::vector<std::array<UM::vec3, 4>> refs_grad_;


	// intern variables 
	std::vector<UM::mat3x3> J_; // per-tet Jacobian matrix = [[JX.x JX.y, JX.z], [JY.x, JY.y, JY.z], [JZ.x, JZ.y, JZ.z]]
	std::vector<UM::mat3x3> K_; // per-tet dual basis: det J = dot J[i] * K[i]
	std::vector<double> det_; // per-tet determinant of the Jacobian matrix
	std::vector<double> vol_; // per-tet determinant of the Jacobian matrix
	double eps_;       // regularization parameter, depends on min(jacobian)

	double detmin_;    // min(jacobian) over all tetrahedra
	int ninverted_; // number of inverted tetrahedra

	double start_eps = 1e-3;

};



//  EXAMPLE OF DEFINITIONS OF VARIABLES FOR TRI AND TET CODE :


class Tri_id_with_lock : public Elliptic_smoother_Variables {
public:
	Tri_id_with_lock(const std::vector<double>& verts, const std::vector<std::array<int, 3>>& triangles, const std::vector<bool>& locks)
	{
		nb_unlocked_ = 0;
		verts_.resize(verts.size());
		locked_.resize(verts.size());

		for (int i = 0; i < (int) verts.size(); i++) verts_[i] = verts[i];
		for (int i = 0; i < (int) verts.size(); i++) locked_[i] = locks[i];

		map_to_vert_.resize(6 * triangles.size());
		reduc_.resize(verts_.size());
		for (int t = 0; t < (int) triangles.size(); t++) {
			for (int tv = 0; tv < 3; tv++) {
				map_to_vert_[6 * t + 2 * tv + 0] = 2 * triangles[t][tv] + 0;
				map_to_vert_[6 * t + 2 * tv + 1] = 2 * triangles[t][tv] + 1;
			}
		}
		for (int i = 0; i < (int) verts_.size(); i++) {
			if (!locked_[i]) {
				reduc_[i] = nb_unlocked_++;
			}
			else
				reduc_[i] = -1;
		}
	}
	inline int nb_of_reducted_variables() { return nb_unlocked_; }
	inline void set_reducted_values(const std::vector<double>& X) {
		for (int i = 0; i < (int) verts_.size(); i++) if (!locked_[i]) {
			verts_[i] = X[reduc_[i]];
		}
	}
	inline void get_reducted_values(std::vector<double>& X) {
		for (int i = 0; i < (int) verts_.size(); i++) if (!locked_[i]) {
			X[reduc_[i]] = verts_[i];
		}
	}

	inline double value(const int i) {
		return verts_[map_to_vert_[i]];
	}
	inline void add2grad(const int i, const double value, std::vector<double>& G) {
		if (!locked_[map_to_vert_[i]]) G[reduc_[map_to_vert_[i]]] += value;
	}
	inline void add2Hessian(const int i, const int j, const double value, std::vector<sparse_term>& H) {
		H.push_back({ {reduc_[map_to_vert_[i]], reduc_[map_to_vert_[j]]}, value });
	}
	void get_verts(std::vector<double>& verts) {
		for (int i = 0; i < (int) verts_.size(); i++) {
			verts[i] = verts_[i];
		}
	}
private:
	int nb_unlocked_;
	std::vector<int> map_to_vert_;
	std::vector<double> verts_;
	std::vector<bool> locked_;
	std::vector<int> reduc_;
};

class Tets_id_with_lock : public Elliptic_smoother_Variables {
public:
	Tets_id_with_lock(const std::vector<double>& verts, const std::vector<std::array<int, 4>>& tets, const std::vector<bool>& locks)
	{
		nb_unlocked_ = 0;
		verts_.resize(verts.size());
		locked_.resize(verts.size());

		for (int i = 0; i < (int) verts.size(); i++) verts_[i] = verts[i];
		for (int i = 0; i < (int) verts.size(); i++) locked_[i] = locks[i];

		map_to_vert_.resize(12 * tets.size());
		reduc_.resize(verts_.size());
		for (int t = 0; t < (int) tets.size(); t++) {
			for (int tv = 0; tv < 4; tv++) {
				for (int d = 0; d < 3; d++) {
					map_to_vert_[12 * t + 3 * tv + d] = 3 * tets[t][tv] + d;
				}
			}
		}
		for (int i = 0; i < (int) verts_.size(); i++) {
			if (!locked_[i]) {
				reduc_[i] = nb_unlocked_++;
			}
			else
				reduc_[i] = -1;
		}
	}
	inline int nb_of_reducted_variables() { return nb_unlocked_; }
	inline void set_reducted_values(const std::vector<double>& X) {
		for (int i = 0; i < (int) verts_.size(); i++) if (!locked_[i]) {
			verts_[i] = X[reduc_[i]];
		}
	}
	inline void get_reducted_values(std::vector<double>& X) {
		for (int i = 0; i < (int) verts_.size(); i++) if (!locked_[i]) {
			X[reduc_[i]] = verts_[i];
		}
	}
	inline double value(const int i) {
		return verts_[map_to_vert_[i]];
	}
	inline void add2grad(const int i, const double value, std::vector<double>& G) {
		if (!locked_[map_to_vert_[i]]) G[reduc_[map_to_vert_[i]]] += value;
	}
	inline void add2Hessian(const int i, const int j, const double value, std::vector<sparse_term>& H) {
		H.push_back({ {reduc_[map_to_vert_[i]], reduc_[map_to_vert_[j]]}, value });
	}
	void get_verts(std::vector<double>& verts) {
		for (int i = 0; i < (int) verts_.size(); i++) {
			verts[i] = verts_[i];
		}
	}
private:
	int nb_unlocked_;
	std::vector<int> map_to_vert_;
	std::vector<double> verts_;
	std::vector<bool> locked_;
	std::vector<int> reduc_;
};