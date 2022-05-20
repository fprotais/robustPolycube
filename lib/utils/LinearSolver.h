#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H
#include <array>
#include <vector>
#include <string>
#include <iostream>
#define DEFAULT_LP_TIMEOUT 300

#include <limits>




extern bool LinearSolver_is_runnable;

bool INITIALISE_LINEARSOLVER(
	const std::string& lpsolve_path,
	const std::string& cplex_path,
	const std::string& glpk_path,
	const std::string& gurobi_cl_path
);



//exemple of use at the end of the file
class LinearSolver {
public:

	constexpr static double infinityplus =  std::numeric_limits<double>::max();
	constexpr static double infinityminus = -std::numeric_limits<double>::max();

	enum optimisation_types : size_t {
		MINIMIZE = 0, MAXIMIZE
	};
	LinearSolver(const optimisation_types& optimisation_type = MINIMIZE, const size_t& solver_time_out = DEFAULT_LP_TIMEOUT);

	void get_actual_solver();



	enum variable_types : size_t {
		FLOAT_VARIABLE = 0, UINT_VARIABLE, BIN_VARIABLE
	};

	enum constraint_types : size_t {
		INEQUALITY_INF = 0, INEQUALITY_SUP, EQUALITY
	};

	size_t add_variable(const double energy_cost, const variable_types variable_type = FLOAT_VARIABLE);
	void set_energy_cost(const size_t variable_id, const double energy_cost);

	size_t add_variables(const std::vector<double>& energy_costs, const std::vector<variable_types>& variables_type = std::vector<variable_types>());

	void fix_variable(const size_t variable_id, const double value = 0);
		
	// default bound is [0, +Inf]
	void add_variable_bound(const size_t variable_id, const double lower_bound, const double upper_bound);

	void begin_constraint(const constraint_types constraint_type = INEQUALITY_INF);
	// adding a variable several time will sum the weights   (TODO : check this)
	void add_constraint_coefficient(const size_t variable_id, const double weight);
	void add_constraint_righthandside(const double value);
	void end_constraint();

	size_t nb_constraints();
	size_t nb_variables();

	const std::string print_problem();
	bool solve();
	const std::string print_solution();
	const std::string print_timings();


	double get_variable(const size_t variable_id);

private:
public:
	optimisation_types _optimisation_type;

	std::vector<double> _variables_value;
	std::vector<double> _variables_energy;
	std::vector<variable_types> _variables_type;
	std::vector<bool> _variables_is_fixed;
	std::vector<size_t> _variables_bounds_index;
	std::vector<std::array<double, 2>> _variables_bounds;

	struct _constraint {
		_constraint(const constraint_types& constraint_type)
			: _ids()
			, _weights()
			, _types()
			, _left_term(0)
			,_type(constraint_type)
		{}
		std::vector<size_t> _ids;
		std::vector<double> _weights;
		std::vector<variable_types> _types;

		double _left_term;
		constraint_types _type;
	};
		
	friend std::ostream& operator<<(std::ostream& os, const _constraint& dt);



	std::vector<_constraint> _constraints;
	bool _is_in_constraint;
	bool _is_solved;
	double _run_time;


	std::string _in_file;
	std::string _out_file;

	size_t _runned_solver;
	size_t _solver_time_out;

	bool run_solver();
	bool read_out_file();


public :
	/*														*
				3 EXAMPLES OF USES OF THIS CLASS						 
		*														*/								 

	inline static void example1() {
		/* Linear Problem
		solving  Min x1 + x2
					st  x1 + 2x2 >= 2,
						x1 >=0,
						x2 >=0,
						x1,x2 in R
			sol : x1 = 0, x2 =1
			https://fr.wikipedia.org/wiki/Optimisation_lin%C3%A9aire */
		std::cout << "================EXAMPLE OF LINEAR PROBLEM SOLVING================" << std::endl;
		std::cout << "Minimise \t x1 + x2 " << std::endl;
		std::cout << "Subject to \t x1 + 2x2 >= 2 " << std::endl;
		std::cout << "\t \t x1  >= 0 " << std::endl;
		std::cout << "\t \t x2  >= 0 " << std::endl;
		std::cout << "And \t\t x1, x2 in R  " << std::endl;

		std::cout << ">>>>>>>>>>> Start of Solver output " << std::endl;
		LinearSolver solver(LinearSolver::MINIMIZE);
		solver.add_variable(1, LinearSolver::FLOAT_VARIABLE); // x1 in R
		solver.add_variable(1, LinearSolver::FLOAT_VARIABLE); // x2 in R

		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(0, 1); // x1
		solver.add_constraint_coefficient(1, 2); // + 2x2
		solver.add_constraint_righthandside(2);  // >= 2
		solver.end_constraint();

		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(0, 1); // x1
		solver.add_constraint_righthandside(0);  // >= 0
		solver.end_constraint();

		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(1, 1); // x2
		solver.add_constraint_righthandside(0);  // >= 0
		solver.end_constraint();

		std::cout << solver.print_problem();
		solver.solve();
		std::cout << solver.print_solution();

		double x1 = solver.get_variable(0);
		double x2 = solver.get_variable(1);

		std::cout << ">>>>>>>>>>> End of Solver output " << std::endl;

		std::cout << "Final solution : " << std::endl;
		std::cout << "x1 = " << x1 << std::endl;
		std::cout << "x2 = " << x2 << std::endl;
		std::cout << "===============================END===============================" << std::endl;

	}

	inline static void example2() {
		/* Mixed Integer Linear Problem
		solving max	y
				st	-x + y <= 1,
					3x + 2y - z <= 8,
					2x + 3y - z <= 8,
					z = 4,
					x,y,z in N
			sol : (1,2) and (2,2)
			https://en.wikipedia.org/wiki/Integer_programming */
		std::cout << "=========EXAMPLE OF MIXED INTEGER LINEAR PROBLEM SOLVING=========" << std::endl;
		std::cout << "Minimise \t y " << std::endl;
		std::cout << "Subject to \t -x + y <= 1 " << std::endl;
		std::cout << "\t \t 3x + 2y - z  <= 8 " << std::endl;
		std::cout << "\t \t 2x + 3y - z  <= 8 " << std::endl;
		std::cout << "\t \t z  = 4 " << std::endl;
		std::cout << "And \t\t x, y, z in N  " << std::endl;

		std::cout << ">>>>>>>>>>> Start of Solver output " << std::endl;

		LinearSolver solver(LinearSolver::MAXIMIZE);
		solver.add_variable(0, LinearSolver::UINT_VARIABLE); // x in N  (>= 0)
		solver.add_variable(1, LinearSolver::UINT_VARIABLE); // y in N  (>= 0)
		solver.add_variable(0, LinearSolver::UINT_VARIABLE); // z in N

		solver.begin_constraint(LinearSolver::INEQUALITY_INF);
		solver.add_constraint_coefficient(0, -1); // -x
		solver.add_constraint_coefficient(1, 1); // y
		solver.add_constraint_righthandside(1);  // <= 1
		solver.end_constraint();

		solver.begin_constraint(LinearSolver::INEQUALITY_INF);
		solver.add_constraint_coefficient(0, 3); // 3x
		solver.add_constraint_coefficient(1, 2); // + 2y
		solver.add_constraint_coefficient(2, -1); // -z
		solver.add_constraint_righthandside(8);  // <= 8
		solver.end_constraint();

		solver.begin_constraint(LinearSolver::INEQUALITY_INF);
		solver.add_constraint_coefficient(0, 2); // 2x
		solver.add_constraint_coefficient(1, 3); // + 3y
		solver.add_constraint_coefficient(2, -1); // -z
		solver.add_constraint_righthandside(8);  // <= 8
		solver.end_constraint();

		solver.fix_variable(2, 4); // z = 4

		std::cout << solver.print_problem();
		solver.solve();
		std::cout << solver.print_solution();

		double x = solver.get_variable(0);
		double y = solver.get_variable(1);
		double z = solver.get_variable(2);

		std::cout << ">>>>>>>>>>> End of Solver output " << std::endl;

		std::cout << "final solution : " << std::endl;
		std::cout << "x = " << x << std::endl;
		std::cout << "y = " << y << std::endl;
		std::cout << "z = " << z << std::endl;

		std::cout << "===============================END===============================" << std::endl;
	}


	inline static void example3() {
		std::cout << "==================EXAMPLE OF ADDING CONSTRAINTS==================" << std::endl;
		std::cout << "We want to find : " << std::endl;
		std::cout << "Minimise \t |x1 - x2| + x2 " << std::endl;
		std::cout << "Subject to \t x1 + x2 >= 2  " << std::endl;
		std::cout << "And \t x1, x2 in N  " << std::endl;
		std::cout << std::endl;

		LinearSolver solver(LinearSolver::MINIMIZE);

		std::cout << "We first try to solve without absolute values: " << std::endl;
		std::cout << "Minimise \t x1   " << std::endl;
		std::cout << "Subject to \t x1 + x2 >= 2  " << std::endl;
		std::cout << "And \t\t x1, x2 in N  " << std::endl;

		std::cout << ">>>>>>>>>>> Start of Solver output " << std::endl;
		solver.add_variables({ {1.,0.} }, { { LinearSolver::UINT_VARIABLE,LinearSolver::UINT_VARIABLE  } }); // x1, x2 \in N

		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(0, 1); // x1
		solver.add_constraint_coefficient(1, 1); // +x2
		solver.add_constraint_righthandside(2); // >= 2
		solver.end_constraint();
		std::cout << solver.print_problem();
		solver.solve();
		std::cout << solver.print_solution();

		double x1 = solver.get_variable(0);
		double x2 = solver.get_variable(1);

		std::cout << ">>>>>>>>>>> End of Solver output " << std::endl;

		std::cout << "This give a first solution : " << std::endl;
		std::cout << "x1 = " << x1 << std::endl;
		std::cout << "x2 = " << x2 << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "Which is obviously not the one we are looking for. Now we linearize the absolute value" << std::endl;
		std::cout << "We first try to solve without absolute values: " << std::endl;
		std::cout << "Minimise \t x3 + x2   " << std::endl;
		std::cout << "Subject to \t x1 + x2 >= 2  " << std::endl;
		std::cout << "\t \t x3 - x1 + x2 >= 0  " << std::endl;
		std::cout << "\t \t x3 + x1 - x2 >= 0  " << std::endl;
		std::cout << "And \t\t x1, x2 in N  " << std::endl;
		std::cout << "\t \t x3 in R  " << std::endl;

		std::cout << ">>>>>>>>>>> Start of Solver output " << std::endl;
		solver.set_energy_cost(0, 0); // 0x1
		solver.set_energy_cost(1, 1); // + 1x2
		solver.add_variables({ 1 }, { LinearSolver::FLOAT_VARIABLE }); // + 1x3

		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(0, -1); // -x1
		solver.add_constraint_coefficient(1, 1); // +x2
		solver.add_constraint_coefficient(2, 1); // +x3
		solver.add_constraint_righthandside(0); // >= 0
		solver.end_constraint();

		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(0, 1); // x1
		solver.add_constraint_coefficient(1, -1); // -x2
		solver.add_constraint_coefficient(2, 1); // +x3
		solver.add_constraint_righthandside(0); // >= 0
		solver.end_constraint();
		std::cout << solver.print_problem();
		solver.solve();
		std::cout << solver.print_solution();

		x1 = solver.get_variable(0);
		x2 = solver.get_variable(1);
		double x3 = solver.get_variable(2);
		std::cout << ">>>>>>>>>>> End of Solver output " << std::endl;

		std::cout << "final solution : " << std::endl;
		std::cout << "x1 = " << x1 << std::endl;
		std::cout << "x2 = " << x2 << std::endl;
		std::cout << "x3 = " << x3 << std::endl;
		std::cout << "Which is all good." << std::endl;

		std::cout << "===============================END===============================" << std::endl;

	}
	inline static void example4() {
		std::cout << "==================EXAMPLE OF BINARY ==================" << std::endl;
		std::cout << "We want to find : " << std::endl;
		std::cout << "Minimise \t x1 + 2 x2 + y1 + 2 y2 " << std::endl;
		std::cout << "Subject to \t bx + by = 1  " << std::endl;
		std::cout << "Subject to \t 2 x1 + x2 = bx  " << std::endl;
		std::cout << "Subject to \t y1 + 2 y2 = by " << std::endl;
		std::cout << "And \t b1, b2 in {0,1}  " << std::endl;
		std::cout << "And \t y1, y2, x1, x2 in [0,+Inf]  " << std::endl;
		std::cout << std::endl;

		LinearSolver solver(LinearSolver::MINIMIZE);


		std::cout << ">>>>>>>>>>> Start of Solver output " << std::endl;
		solver.add_variable(0, LinearSolver::BIN_VARIABLE); // b1
		solver.add_variable(0, LinearSolver::BIN_VARIABLE); // b2

		solver.add_variable(1, LinearSolver::FLOAT_VARIABLE); // x1
		solver.add_variable(2, LinearSolver::FLOAT_VARIABLE); // x2

		solver.add_variable(1, LinearSolver::FLOAT_VARIABLE); // y1
		solver.add_variable(2, LinearSolver::FLOAT_VARIABLE); // y2

		solver.begin_constraint(LinearSolver::EQUALITY);
		solver.add_constraint_coefficient(0, 1); 
		solver.add_constraint_coefficient(1, 1);
		solver.add_constraint_righthandside(1); 
		solver.end_constraint();

		solver.begin_constraint(LinearSolver::EQUALITY);
		solver.add_constraint_coefficient(2, 1);
		solver.add_constraint_coefficient(3, 2);
		solver.add_constraint_coefficient(0, -1); 
		solver.end_constraint();
		solver.begin_constraint(LinearSolver::EQUALITY);
		solver.add_constraint_coefficient(4, 2);
		solver.add_constraint_coefficient(5, 1);
		solver.add_constraint_coefficient(1, -1); 
		solver.end_constraint();


		std::cout << solver.print_problem();
		solver.solve();
		std::cout << solver.print_solution();


		double bx = solver.get_variable(0);
		double by = solver.get_variable(1);
		double x1 = solver.get_variable(2);
		double x2 = solver.get_variable(3);
		double y1 = solver.get_variable(4);
		double y2 = solver.get_variable(5);
		std::cout << ">>>>>>>>>>> End of Solver output " << std::endl;

		std::cout << "final solution : " << std::endl;
		std::cout << "bx = " << bx << std::endl;
		std::cout << "by = " << by << std::endl;
		std::cout << "x1 = " << x1 << std::endl;
		std::cout << "x2 = " << x2 << std::endl;
		std::cout << "y1 = " << y1 << std::endl;
		std::cout << "y2 = " << y2 << std::endl;

		std::cout << "===============================END===============================" << std::endl;

	}
};


#endif // !LINEAR_SOLVER_H
