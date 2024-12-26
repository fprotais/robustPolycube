#include "elliptic_smoothing.h"

#include <ultimaille/syntactic-sugar/HLBFGS_wrapper.h>

#include <cmath>
#include <iostream>
#include "algorithm"

#define FOR(i, n) for(int i = 0; i < n; i++)

#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = std::vector<double>(omp_orig.size(), 0))
#endif


inline double chi(double eps, double det) {
	const double eps2 = eps * eps;
    if (det > 0)
        return (det + std::sqrt(eps2 + det * det)) * .5;
    return .5 * eps2 / (std::sqrt(eps2 + det * det) - det);
}

inline double chi_deriv(double eps, double det) {
    return .5 + det / (2. * std::sqrt(eps * eps + det * det));
}

inline void display_options(const smoother_options& opt) {
    std::cerr << "Options: " << std::endl;
    std::cerr << "-theta = " << opt.theta << std::endl;
    std::cerr << "-maxiter = " << opt.maxiter << std::endl;
    std::cerr << "-bfgs_threshold = " << opt.bfgs_threshold << std::endl;
    std::cerr << "-bfgs_maxiter = " << opt.bfgs_maxiter << std::endl;
    std::cerr << "-debug = " << opt.debug << std::endl;
    std::cerr << "-eps_from_theorem = " << opt.eps_from_theorem << std::endl;
    std::cerr << "-stopping_when_static = " << opt.stopping_when_static << std::endl;
    std::cerr << "-static_threshold = " << opt.static_threshold << std::endl;
}



//-----------------------------------------------------------------------------------//
//                                                                                   //
//                                        2D                                         //
//                                                                                   //
//-----------------------------------------------------------------------------------//


void Elliptic_smoother_2D::set_triangle_ref(const int t, const std::array<double, 6> &ref) {
    assert(t < N_tri_);
    mat22 M = { ref[2] - ref[0], ref[4] - ref[0], ref[3] - ref[1], ref[5] - ref[1] };
    double detM = (M[0] * M[3] - M[1] * M[2]);
    area_[t] = 0.5 * detM;
    double invdetM = 1. / detM;
    mat22 invM = { invdetM * M[3], -invdetM * M[1], -invdetM * M[2], invdetM * M[0] };
    refs_grad_[t] = { -invM[0] - invM[2], -invM[1] - invM[3], invM[0], invM[1], invM[2], invM[3] };
}
void Elliptic_smoother_2D::set_triangle_ref(const std::vector<std::array<double, 6>> &refs) {
    assert(refs.size() == N_tri_);
    FOR(t, N_tri_) set_triangle_ref(t, refs[t]);
}
void Elliptic_smoother_2D::set_area(const int t, const double area) {
    assert(t < N_tri_);
    area_[t] = area;
}

void Elliptic_smoother_2D::init() {
    refs_grad_.resize(N_tri_);
    area_.resize(N_tri_);
    FOR(t, N_tri_) {
        refs_grad_[t] = { { 0,-1, std::sqrt(3.) / 2.,.5, -std::sqrt(3.) / 2.,.5 } };
        area_[t] = std::sqrt(3) / 4;
    }
    J_.resize(N_tri_);
    K_.resize(N_tri_);
    det_.resize(N_tri_);
}

void Elliptic_smoother_2D::evaluate_jacobian() {
    if (options_.debug > 3) std::cerr << "evaluate the jacobian...";
    detmin_ = std::numeric_limits<double>::max();
    ninverted_ = 0;
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for reduction(+:ninverted_) reduction(min:detmin_)
#endif
    FOR(t, N_tri_) {
        J_[t] = { 0,0,0,0 };
        FOR(i, 3) FOR(d, 2) {
			const double v = var_.value(6*t+2*i+d);
			FOR(dJ, 2) J_[t][2*d + dJ] += refs_grad_[t][2*i+dJ] * v;
		}

        det_[t] = J_[t][0] * J_[t][3] - J_[t][1] * J_[t][2];
		if(det_[t] < detmin_) detmin_ = det_[t];
		if(det_[t] <= 0.) ++ ninverted_;

        K_[t] = { { J_[t][3], -J_[t][2], -J_[t][1], J_[t][0] } };// dual basis
    }
    if (options_.debug > 3) std::cerr << "ok" << std::endl;
}

double Elliptic_smoother_2D::evaluate_energy() {
    evaluate_jacobian();
    double E = 0;
    FOR(t, N_tri_) {
        double chi_ = chi(eps_, det_[t]);
        double f = 0;
        FOR(i, 4) f += J_[t][i] * J_[t][i];
        double g = 1 + det_[t] * det_[t];
        E += ((1. - options_.theta) * f + options_.theta * g) * area_[t] / chi_;
    }
    return E;
}

void Elliptic_smoother_2D::standard_elliptic_energy(const std::vector<double>& X, double& F, std::vector<double>& G) {
    evaluate_jacobian();

#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for reduction(+:F) reduction(vec_double_plus:G)
#endif
    FOR(t, N_tri_) {
        double c1 = chi(eps_, det_[t]);
        double c2 = chi_deriv(eps_, det_[t]);

        double f = 0;
        FOR(i, 4) f += J_[t][i] * J_[t][i];
        f /= c1;
        double g = (1 + det_[t] * det_[t]) / c1;
        F += ((1. - options_.theta) * f + options_.theta * g) * area_[t];

        const double df_mul_a = 2. / c1;
        const double df_mul_b = f * c2 / c1;
        const double dg_mul = (2 * det_[t] - g * c2) / c1;

        FOR(d, 2) {
            double a1 = J_[t][2 * d], a2 = J_[t][2 * d + 1]; // tangent basis
            double b1 = K_[t][2 * d], b2 = K_[t][2 * d + 1]; // dual basis
            double dfda1 = a1 * df_mul_a - b1 * df_mul_b, dfda2 = a2 * df_mul_a - b2 * df_mul_b;
            double dgda1 = b1 * dg_mul, dgda2 = b2 * dg_mul;
            FOR(i, 3) {
                double gradi = (dfda1 * (1. - options_.theta) + dgda1 * options_.theta) * refs_grad_[t][2 * i] 
                             + (dfda2 * (1. - options_.theta) + dgda2 * options_.theta) * refs_grad_[t][2 * i + 1];
                gradi *= area_[t];
                var_.add2grad(6 * t + 2 * i + d, gradi, G);
            }
        }
    }
}

void Elliptic_smoother_2D::standard_elliptic_energy_w_hessian(const std::vector<double>& X, double& F, std::vector<double>& G, std::vector<Elliptic_smoother_Variables::sparse_term>& H) {
	standard_elliptic_energy(X, F, G);
}

bool Elliptic_smoother_2D::run_lbfgs(std::vector<double>& X) {
	bool first = true;
	double F0, F1;
    const STLBFGS::Optimizer::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
        std::fill(G.begin(), G.end(), 0);
        F = 0;
        var_.set_reducted_values(X);
        energy(X, F, G);
		if(first) {
			F0 = F1 = F;
			first = false;
		} else F1 = std::min(F1, F);
    };
    STLBFGS::Optimizer opt = { func };
    //opt.gtol = options_.bfgs_threshold;
    opt.ftol = options_.bfgs_threshold;
    //opt.maxiter = options_.bfgs_maxiter;
    //opt.verbose = options_.debug;;
    opt.invH.history_depth = 20;
    //UM::LBFGS_Optimizer opt(func);

    //if(!options_.eps_from_theorem && detmin_ > 0.) {
    //    opt.gtol = std::min(options_.bfgs_threshold, .1 * options_.static_threshold);
    //    opt.maxiter = options_.bfgs_maxiter * std::max(1, options_.maxiter / 2);
    //} else {
        opt.gtol = options_.bfgs_threshold;
        opt.maxiter = options_.bfgs_maxiter;
    //}
    opt.verbose = options_.debug > 0;


    opt.run(X);
    if (options_.debug > 0) std::cerr << " F : " << F0 << "  --->  " << F1 << "    ||   ";
    return true;// (F0 - F1) / F1 > 2.e-7;
}

bool Elliptic_smoother_2D::run_newton(std::vector<double>& X) {
    return true;
}

bool Elliptic_smoother_2D::go() {
    if (options_.debug > 0) {
        std::cerr << "==== Running Elliptic smoother 2d. ====" << std::endl;
        display_options(options_);
    }
    if (var_.nb_of_reducted_variables() == 0) {
        std::cerr << "No variables to optimize" << std::endl;
        return true;
    }
    double total_area_ = 0;
    FOR(i, N_tri_) total_area_ += area_[i];
    FOR(i, N_tri_) area_[i] /= total_area_;

    evaluate_jacobian();
    const double e0 = 1e-3;
    if (options_.eps_from_theorem)
        eps_ = start_eps;
    else if (options_.barrier)
        eps_ = 1e-7;
    else
        eps_ = detmin_ > 0 ? .5*e0 : std::sqrt(e0*e0 + 0.04 * detmin_ * detmin_);
    std::vector<double> X(var_.nb_of_reducted_variables());
    var_.get_reducted_values(X);
    bool nullStep = false;

    FOR(iter, options_.maxiter) {
        const double E_prev = evaluate_energy();
        const double detmin_prev = detmin_;
        if (options_.debug > 0) std::cerr << "iteration #" << iter << ":    eps: " << eps_ << " detmin: " << detmin_ << " ninv: " << ninverted_ << std::endl;

		bool better;
        if (options_.use_newton) better = run_newton(X);
        else better = run_lbfgs(X);
        var_.set_reducted_values(X);

        const double E = evaluate_energy();
        if (options_.debug > 0) std::cerr << " E : " << E_prev << "  --->  " << E << std::endl;
        if (options_.eps_from_theorem) {
            const double sigma = std::max(1. - E / E_prev, 1e-1);
            if (detmin_ >= 0) eps_ *= (1 - sigma);
            else {
				const double det_eps_norm = std::sqrt(detmin_*detmin_ + eps_*eps_);
                eps_ *= 1 - (sigma * det_eps_norm) / (std::abs(detmin_) + det_eps_norm);
			}
        } else if (detmin_prev > 0. && detmin_ > 0.) {
            std::cerr << "Stopping as detmin > 0 while not using the eps from theorem" << std::endl;
            break;
        } else if(!options_.barrier) {
            eps_ = std::min(.995*eps_, detmin_ > 0 ? .5*e0 : std::sqrt(e0*e0 + 0.04 * detmin_*detmin_));
        }
        if (detmin_ > 0 || options_.stopping_when_static) {
            if ((E_prev - E) / E < options_.static_threshold) break;
        } else if (!better) {
            if(nullStep) break;
            nullStep = true;
        } else nullStep = false;
    }
    if (options_.debug > 0) std::cerr << "E: " << evaluate_energy() << " detmin: " << detmin_ << " ninv: " << ninverted_ << std::endl;
    return !ninverted_;
}



//-----------------------------------------------------------------------------------//
//                                                                                   //
//                                        3D                                         //
//                                                                                   //
//-----------------------------------------------------------------------------------//


void Elliptic_smoother_3D::set_tet_ref(const int t, const std::array < UM::vec3, 4> tet_ref) {
    assert(t < N_tets_);
    UM::mat3x3 M = { tet_ref[1] - tet_ref[0], tet_ref[2] - tet_ref[0], tet_ref[3] - tet_ref[0] };
    UM::mat3x3 invM = M.invert();
    double detM = M.det();
    vol_[t] = 1./6. * detM;
    invM = invM.transpose();
    refs_grad_[t] = { -invM[0] - invM[1] - invM[2], invM[0], invM[1], invM[2] };
}
void Elliptic_smoother_3D::set_tet_ref(const std::vector<std::array < UM::vec3, 4> > tets_ref) {
    assert(tets_ref.size() == N_tets_);
    for (int t = 0; t < N_tets_; t++) set_tet_ref(t, tets_ref[t]);
}

void Elliptic_smoother_3D::init() {
    refs_grad_.resize(N_tets_);
    vol_.resize(N_tets_);
    FOR(t, N_tets_) {
        constexpr double a = 0.70710678118; // invsqrt of 2, not possible in constexpr...
        // ref tet is : [0,0,0], [a,0,a], [a,a,0], [0,a,a]
        refs_grad_[t] = { UM::vec3(-a,-a,-a), UM::vec3(a,a,-a), UM::vec3(-a,a,a), UM::vec3(a,-a,a) };
        vol_[t] = 1./6.*a;
    }
    J_.resize(N_tets_);
    K_.resize(N_tets_);
    det_.resize(N_tets_);
}
inline UM::mat3x3 dual_basis(const UM::mat3x3& J) {
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

void Elliptic_smoother_3D::evaluate_jacobian() {
    if (options_.debug > 3) std::cerr << "evaluate the jacobian...";
    detmin_ = std::numeric_limits<double>::max();
    ninverted_ = 0;
    FOR(t, N_tets_) {
        J_[t] = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
        FOR(i, 4) FOR(d, 3)
            J_[t][d] += refs_grad_[t][i] * var_.value(12 * t + 3 * i + d);

        det_[t] = J_[t].det();
        detmin_ = std::min(detmin_, det_[t]);
        ninverted_ += (det_[t] <= 0);

        K_[t] = dual_basis(J_[t]);
    }
    if (options_.debug > 3) std::cerr << "ok" << std::endl;

}
double Elliptic_smoother_3D::evaluate_energy() {
    evaluate_jacobian();
    double E = 0;
    FOR(t, N_tets_) {
        double c = chi(eps_, det_[t]);
        double f = (J_[t][0] * J_[t][0] + J_[t][1] * J_[t][1] + J_[t][2] * J_[t][2]) / std::pow(c, 2. / 3.);
        double g = (1 + det_[t] * det_[t]) / c;
        E += ((1 - options_.theta) * f + options_.theta * g) * vol_[t];
    }
    return E;
}

void Elliptic_smoother_3D::standard_elliptic_energy(const std::vector<double>&, double& F, std::vector<double>& G) {
    F += evaluate_energy();
    FOR(t, N_tets_) {
        double c1 = chi(eps_, det_[t]);
        double c2 = std::pow(c1, 2. / 3.);
        double c3 = chi_deriv(eps_, det_[t]);
        double f = (J_[t][0] * J_[t][0] + J_[t][1] * J_[t][1] + J_[t][2] * J_[t][2]) / c2;
        double g = (1 + det_[t] * det_[t]) / c1;


        FOR(d, 3) {
            UM::vec3 dfda = J_[t][d] * (2. / c2) - K_[t][d] * ((2. * f * c3) / (3. * c1));
            UM::vec3 dgda = K_[t][d] * ((2 * det_[t] - g * c3) / c1);

            FOR(tc, 4) {
                double gradi = (dfda * (1. - options_.theta) + dgda * options_.theta) * refs_grad_[t][tc] * vol_[t];
                var_.add2grad(12 * t + 3 * tc + d, gradi, G);
            }
        }
    }

}

void Elliptic_smoother_3D::standard_elliptic_energy_w_hessian(const std::vector<double>& X, double& F, std::vector<double>& G, std::vector<Elliptic_smoother_Variables::sparse_term>& H) {
    F += evaluate_energy();
    FOR(t, N_tets_) {
        double c1 = chi(eps_, det_[t]);
        double c2 = std::pow(c1, 2. / 3.);
        double c3 = chi_deriv(eps_, det_[t]);
        double f = (J_[t][0] * J_[t][0] + J_[t][1] * J_[t][1] + J_[t][2] * J_[t][2]) / c2;
        double g = (1 + det_[t] * det_[t]) / c1;


        FOR(d, 3) {
            UM::vec3 dfda = J_[t][d] * (2. / c2) - K_[t][d] * ((2. * f * c3) / (3. * c1));
            UM::vec3 dgda = K_[t][d] * ((2 * det_[t] - g * c3) / c1);

            FOR(tc, 4) {
                double gradi = (dfda * (1. - options_.theta) + dgda * options_.theta) * refs_grad_[t][tc] * vol_[t];
                var_.add2grad(12 * t + 3 * tc + d, gradi, G);
            }
        }
    }

}

bool Elliptic_smoother_3D::run_lbfgs(std::vector<double>& X) {
    const STLBFGS::Optimizer::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
        std::fill(G.begin(), G.end(), 0);
        F = 0;

        var_.set_reducted_values(X);
        energy(X, F, G);
    };
    STLBFGS::Optimizer opt = { func };

    opt.gtol = options_.bfgs_threshold;
    opt.ftol = options_.bfgs_threshold;
    opt.maxiter = options_.bfgs_maxiter;
    opt.verbose = options_.debug;
    opt.invH.history_depth = 5;
    //UM::LBFGS_Optimizer opt(func);
    //opt.gtol = options_.bfgs_threshold;
    //opt.maxiter = options_.bfgs_maxiter;
    //opt.verbose = options_.debug;
    opt.run(X);
    return true;
}

bool Elliptic_smoother_3D::run_newton(std::vector<double>& X) {
    return true;
}

bool Elliptic_smoother_3D::go() {
    if (options_.debug > 0) {
        std::cerr << "==== Running Elliptic smoother 3d. ====" << "\n";
        display_options(options_);
    }
    if (var_.nb_of_reducted_variables() == 0) {
        std::cerr << "No variables to optimize" << "\n";
        return true;
    }
    double total_vol_ = 0; // so that energy is always roughly the same
    FOR(i, N_tets_) total_vol_ += vol_[i];
    if (options_.debug > 0) std::cerr << "total vol = " << total_vol_ << "\n";
    FOR(i, N_tets_) vol_[i] = vol_[i] / total_vol_;
    evaluate_jacobian();
    double e0 = 1e-3;
    if (options_.eps_from_theorem) {
        eps_ = start_eps;
        double sigma = 0.5;
        double mu = (1 - sigma) * chi(eps_, detmin_);
        if (detmin_ < mu)
            eps_ = 2 * std::sqrt(mu * (mu - detmin_));
        else eps_ = 1e-10;
    }
    if (options_.barrier) {
        e0 = 1e-7;
        eps_ = 1e-7;
    }
    FOR(iter, options_.maxiter) {
        if (options_.debug > 0) std::cerr << "iteration #" << iter << "\n";
        if (!options_.eps_from_theorem) {
            if (iter && iter % 10 == 0 && e0 > 1e-10) e0 /= 2.;
            eps_ = detmin_ > 0 ? e0 : std::sqrt(e0 * e0 + 0.04 * detmin_ * detmin_);
        }
        if (options_.debug > 0) std::cerr << "E: " << evaluate_energy() << " eps: " << eps_ << " detmin: " << detmin_ << " ninv: " << ninverted_ << std::endl;


        std::vector<double> X(var_.nb_of_reducted_variables());
        var_.get_reducted_values(X);
        double E_prev = evaluate_energy();
        if (options_.use_newton) run_newton(X);
        else run_lbfgs(X);
        var_.set_reducted_values(X);

        double E = evaluate_energy();
        if (options_.debug > 0) std::cerr << " E : " << E_prev << "  --->  " << E << "\n";
        if (options_.eps_from_theorem) {
            double sigma = std::max(1. - E / E_prev, 0.1);
            double mu = (1 - sigma) * chi(eps_, detmin_);
            if (detmin_ < mu)
                eps_ = 2 * std::sqrt(mu * (mu - detmin_));
            else eps_ = 1e-10;
        }
        iter_call_back(iter);
        if (detmin_ > 0 && std::abs(E_prev - E) / E < options_.static_threshold) break;
        if (options_.stopping_when_static && std::abs(E_prev - E) / E < options_.static_threshold) {
            std::cerr << "Stopping as Energy is static with threshold: " << options_.static_threshold << std::endl;
            break;
        }

    }
    if (options_.debug > 0) std::cerr << "E: " << evaluate_energy() << " detmin: " << detmin_ << " ninv: " << ninverted_ << std::endl;
    return !ninverted_;
}
