#include "LinearSolver.h"

#include <ctime>
#include <sstream>
#include <iostream>
#include <fstream>

#include <filesystem>
//geogram dependancy

const size_t NOT_AN_ID = size_t(-1);


#define TALK(x) std::cerr << x << std::endl

inline std::string disp_file(const char* file_in, int line) {
	std::string file(file_in);
	return std::filesystem::path(file).stem().string() + std::filesystem::path(file).extension().string() + ":" + std::to_string(line);
}
#define UNWANTED_BEHAVIOR(x) {TALK(x);TALK(disp_file(__FILE__, __LINE__)); exit(1);}

static std::string lp_solve_binary_path = "";
static std::string cplex_binary_path = "";
static std::string glpk_binary_path = "";
static std::string gurobi_cl_binary_path = "";

bool LinearSolver_is_runnable = false;

bool INITIALISE_LINEARSOLVER(
	const std::string& lpsolve_path,
	const std::string& cplex_path,
	const std::string& glpk_path,
	const std::string& gurobi_cl_path
) {
	std::string file_name = "tmpfile.txt";
#ifdef WIN32
	std::string redirect_tmp_file = " > " + file_name + " 2>&1";
#else
	std::string redirect_tmp_file = ">  " + file_name;
#endif // win32
	lp_solve_binary_path = lpsolve_path;
	cplex_binary_path = cplex_path;
	glpk_binary_path = glpk_path;
	gurobi_cl_binary_path = gurobi_cl_path;

	{
		if (!std::system(std::string(lpsolve_path + " -h " + redirect_tmp_file).c_str())) {
			TALK("--found lp_solve : " << lpsolve_path);
			LinearSolver_is_runnable = true;
		}
	}
	{
		if (!std::system(std::string(glpk_path + " --version " + redirect_tmp_file).c_str())) {
			TALK("--found glpk : " << glpk_path);
			LinearSolver_is_runnable = true;
		}
	}
	{
		if (!std::system(std::string(cplex_path + " -c \"help\" " + redirect_tmp_file).c_str())) {
			TALK("--found cplex : " << cplex_path);
			LinearSolver_is_runnable = true;
		}
		std::filesystem::remove("cplex.log");
	}
	{
		if (!std::system(std::string(gurobi_cl_path + " -h " + redirect_tmp_file).c_str())) {
			TALK("--found gurobi : " << gurobi_cl_path);
			LinearSolver_is_runnable = true;
		}
		std::filesystem::remove("gurobi.log");
	}
	std::filesystem::remove(file_name);
	return LinearSolver_is_runnable;
}





inline static bool is_int(const double value) {
	return (value == std::round(value));
}

inline static bool is_uint(const double value) {
	return is_int(value);
}
inline static bool is_BIN(const double value) {
	return (value == 0 || value == 1);
}

inline std::ostream& operator<<(std::ostream& os, const LinearSolver::_constraint& cons)
{
	for(size_t v = 0; v< cons._ids.size(); ++v){
		if (cons._weights[v] <= 0)
			os << " - " << std::abs(cons._weights[v]);
		else 
			os << " + " << cons._weights[v];
		if (cons._types[v] == LinearSolver::FLOAT_VARIABLE)
			os << " F_" << cons._ids[v];
		else if (cons._types[v] == LinearSolver::UINT_VARIABLE)
			os << " UI_" << cons._ids[v];
		else if (cons._types[v] == LinearSolver::BIN_VARIABLE)
			os << " B_" << cons._ids[v];
		else
			UNWANTED_BEHAVIOR("Unknown variable type, this is not supposed to happens, please investigate.");
	}
	if (cons._type == LinearSolver::INEQUALITY_INF)
		os << " <= ";
	else if (cons._type == LinearSolver::INEQUALITY_SUP)
		os << " >= ";		
	else 
		os << " = ";
	os << cons._left_term;
	return os;
}

LinearSolver::LinearSolver(const optimisation_types& optimisation_type, const size_t& solver_time_out)
	: _optimisation_type(optimisation_type)
	, _is_in_constraint(false)
	, _is_solved(false)
	, _in_file("polycube.lp")
	, _out_file("polycube.sol")
	, _solver_time_out(solver_time_out)
{
	TALK("Initialising LinearSolver :");
	get_actual_solver();
}

void LinearSolver::get_actual_solver() {
	std::string file_name =  "tmpfile.txt";
#ifdef WIN32
	std::string redirect_tmp_file = " > " +   file_name + " 2>&1";
#else
	std::string redirect_tmp_file = ">  " +  file_name;
#endif // win32


	{
		if (!std::system(std::string(lp_solve_binary_path + " -h " + redirect_tmp_file).c_str())) {
			TALK("--found lp_solve : " << lp_solve_binary_path);
			_runned_solver = 1;
		}
	}
	{
		if (!std::system(std::string(glpk_binary_path + " --version " + redirect_tmp_file).c_str())) {
			TALK("--found glpk : " << glpk_binary_path);
			_runned_solver = 3;
		}
	}
	{
		if (!std::system(std::string(cplex_binary_path + " -c \"help\" " + redirect_tmp_file).c_str())) {
			TALK("--found cplex : "<< cplex_binary_path);
			_runned_solver = 2;
		}
		std::filesystem::remove("cplex.log");

	}
	{
		if (!std::system(std::string(gurobi_cl_binary_path + " -h " + redirect_tmp_file).c_str())) {
			TALK("--found gurobi : " << gurobi_cl_binary_path);
			_runned_solver = 4;
		}
		std::filesystem::remove("gurobi.log");
	}
	std::filesystem::remove(file_name);
	TALK("chosing best solver : ");
	switch (_runned_solver)
	{
	case 1:
		TALK( " ======= Using LP_SOLVE  ======= \n path is : " + lp_solve_binary_path);
		break;
	case 2:
		TALK( "======= Using CPLEX  =======  \n path is : " + cplex_binary_path);
		break;
	case 3:
		TALK("======= Using GLPK - GLPSOL  ======= \n path is : " + glpk_binary_path);
		break;
	case 4:
		TALK("======= Using GUROBI_CL  ======= \n path is : " + gurobi_cl_binary_path);
		break;
	default:
		UNWANTED_BEHAVIOR("There was no solver detected, one must be initialised before running. \n"
			"To do so, you must install a solver and set the variable in hexdom. You can set the variables in ~/.hexdom_default.command, in your command file, or in hexdom_config_directory.cpp\n"
			"The compatible solvers are CPLEX, LP_SOLVE and GLKP-GLPSOL.\n"
			"- CPLEX : set cplex_binary=cplex \n"
			"- GUROBI : set gurobi_cl_binary=gurobi_cl \n"
			"- GLPK-GPLSOL : set glpk_binary=$PATH/glpk/glpsol.exe \n"
			"- LP_SOLVE : set lp_solve_binary=$PATH/lp_solve/lp_solve.exe \n"
			"For UNIX : the install_everything.sh script will automaticly install glpk for you.\n"
			"For Windows : if you decompress softs.zip into hexdom, there will be lp_solve included\n"
			"We chose what we think is the best available solver for the situation. (GUROBI > CPLEX > all)");
	}
}

size_t LinearSolver::add_variable(const double energy_cost, const variable_types variable_type) {
	_is_solved = false;
	_variables_value.push_back(0);
	_variables_energy.push_back(energy_cost);
	_variables_type.push_back(variable_type);
	_variables_is_fixed.push_back(false);
	_variables_bounds_index.push_back(NOT_AN_ID);
	return _variables_type.size() - 1;
}
void LinearSolver::set_energy_cost(const size_t variable_id, const double energy_cost) {
	if (variable_id > nb_variables())
		UNWANTED_BEHAVIOR("out of bound of variable  : " << variable_id << " > " << nb_variables() - 1);
	_variables_energy[variable_id] = energy_cost;
}

size_t LinearSolver::add_variables(const std::vector<double>& energy_costs, const std::vector<variable_types>& variables_type) {
	std::vector<variable_types> types = variables_type;
	if (types.size() == 0)
		types = std::vector<variable_types>(energy_costs.size(), FLOAT_VARIABLE);
	if (types.size() != energy_costs.size()) UNWANTED_BEHAVIOR("You tried to add " << energy_costs.size() << " variables with " << types.size() << " types. ");
	_is_solved = false;
	size_t nb_of_element = nb_variables();

	_variables_value.resize(nb_of_element + energy_costs.size(), 0);
	_variables_bounds_index.resize(nb_of_element + energy_costs.size(), NOT_AN_ID);
	_variables_energy.resize(nb_of_element + energy_costs.size());
	_variables_type.resize(nb_of_element + energy_costs.size());
	for (size_t i = 0; i < energy_costs.size(); i++) {
		_variables_energy[nb_of_element + i] = energy_costs[i];
		_variables_type[nb_of_element + i] = variables_type[i];
	}
	_variables_is_fixed.resize(nb_of_element + energy_costs.size(), false);

	return _variables_type.size() - 1;
}




void LinearSolver::fix_variable(const size_t variable_id, const double value ) {
	if (variable_id > nb_variables()) 
		UNWANTED_BEHAVIOR("out of bound of variable  : " << variable_id << " > " << nb_variables() - 1);
	_is_solved = false;
	if (_variables_type[variable_id] != FLOAT_VARIABLE && value < 0)
		UNWANTED_BEHAVIOR("You tried to fixed the variable " << variable_id << " which is BIN/unsigned int to a signed value\n This is unsolvable, aborting. ")

	if (_variables_is_fixed[variable_id] && value != _variables_value[variable_id]) 
		UNWANTED_BEHAVIOR("You tried to fixed the variable " << variable_id << " twice to different value \n This is not possible in Linear Programming, aborting. ")
	else {
		_variables_is_fixed[variable_id] = true;

		if (_variables_type[variable_id] == UINT_VARIABLE && !is_uint(value)) {
			TALK("You initialised a unsigned int variable to a float value, rounding ...");
			_variables_value[variable_id] = std::round(value);
		}
		else if (_variables_type[variable_id] == BIN_VARIABLE && !is_BIN(value)) {
			TALK("You initialised a BIN variable to a float value, rounding ...");
			_variables_value[variable_id] = (double)(bool)std::round(value);
		}
		else 
			_variables_value[variable_id] = value;
	}

}
void LinearSolver::add_variable_bound(const size_t variable_id, const double lower_bound, const double upper_bound) {
	if (lower_bound > upper_bound) UNWANTED_BEHAVIOR("You tried to add an uncoherent bound for variable" << variable_id << ", implying : " << lower_bound << " > " << upper_bound << "\n This is not possible in Linear Programming, aborting.");
	if (_variables_type[variable_id] == BIN_VARIABLE) UNWANTED_BEHAVIOR("Bounds for binaries variables are bad practices. Please use the fix_variable feature.\n Aborting.");
	if (_variables_bounds_index[variable_id] == NOT_AN_ID) {
		_variables_bounds_index[variable_id] = _variables_bounds.size();
		_variables_bounds.push_back({ lower_bound, upper_bound });
	}
	else 
		_variables_bounds[_variables_bounds_index[variable_id]] = { lower_bound, upper_bound };

}

void LinearSolver::begin_constraint(const constraint_types constraint_type ) {
	_is_solved = false;
	if (_is_in_constraint) UNWANTED_BEHAVIOR("You tried to begin a constraint when in a constraint");
	_is_in_constraint = true;
	_constraints.push_back(_constraint(constraint_type));
}
void LinearSolver::add_constraint_coefficient(const size_t variable_id, const double weight) {
	if (!_is_in_constraint) UNWANTED_BEHAVIOR("You tried to modify a constraint when out of a constraint");
	if (variable_id > nb_variables())
		UNWANTED_BEHAVIOR("out of bound of variable  : " << variable_id << " > " << nb_variables() - 1);
	_constraints.back()._ids.push_back(variable_id);
	_constraints.back()._weights.push_back(weight);
	_constraints.back()._types.push_back(_variables_type[variable_id]);

}
void LinearSolver::add_constraint_righthandside(const double value) {
	if (!_is_in_constraint) UNWANTED_BEHAVIOR("You tried to modify a constraint when out of a constraint");
	_constraints.back()._left_term = value;
}
void LinearSolver::end_constraint() {
	if (!_is_in_constraint) UNWANTED_BEHAVIOR("You tried to end a constraint when out of a constraint");

	// some weak checking for warnings
	bool all_int = true;
	for (size_t v = 0; v < _constraints.back()._ids.size(); ++v) 
		all_int = all_int && (_constraints.back()._types[v] == UINT_VARIABLE || _constraints.back()._types[v] == BIN_VARIABLE);
	if (all_int && _constraints.back()._type == EQUALITY && !is_int(_constraints.back()._left_term) ) {
		bool all_weight_int = true;
		for (size_t v = 0; v < _constraints.back()._weights.size(); ++v) 
			all_weight_int = all_weight_int && is_uint(_constraints.back()._weights[v]);
		if (!all_weight_int || (all_weight_int && is_uint(_constraints.back()._left_term)))
			TALK("The constraint looks ill posed, be careful with integers and equalities : \n "<< _constraints.back()) << std::endl;
	}

	_is_in_constraint = false;
}

size_t LinearSolver::nb_constraints() {
	return (size_t)_constraints.size();
}
size_t LinearSolver::nb_variables() {
	return _variables_value.size();
}

inline size_t sizeOf_sst(std::stringstream& oss) {
	oss.seekg(0, std::ios::end);
	return oss.tellg();
}


const std::string LinearSolver::print_problem() {
	std::ostringstream  problem;
	if (_optimisation_type == MINIMIZE) {
		problem << "Minimize" << std::endl;;
	}
	else if (_optimisation_type == MAXIMIZE) {
		problem << "Maximize" << std::endl;;
	}
	problem << " obj : ";
	std::stringstream size_ensurer;
	for (size_t v = 0; v < nb_variables(); ++v){
 		if (_variables_energy[v] < 0)
			size_ensurer << " - " << std::abs(_variables_energy[v]);
		else
			size_ensurer << " + " << _variables_energy[v];
		if (_variables_type[v] == LinearSolver::FLOAT_VARIABLE)
			size_ensurer << " F_" << v;
		else if (_variables_type[v] == LinearSolver::UINT_VARIABLE)
			size_ensurer << " UI_" << v;
		else if (_variables_type[v] == LinearSolver::BIN_VARIABLE)
			size_ensurer << " B_" << v;
		else 
			UNWANTED_BEHAVIOR("Unknown variable type, this is not supposed to happens, please investigate.");
		if (sizeOf_sst(size_ensurer) > 5000) {
			//std::cerr << "Push : " << sizeOf_sst(size_ensurer)  << " -> " << problem.str().size() << std::endl;

			problem << size_ensurer.str() << std::endl;
			size_ensurer.str(std::string());

		}
	}
	problem << size_ensurer.str();
	problem << std::endl;
	problem << "Subject To " << std::endl;
	for (size_t c = 0; c < nb_constraints(); ++c) {
		problem << _constraints[c] << std::endl;
	}

	problem << "\\ fixed variables " << std::endl;
	for (size_t v = 0; v < nb_variables(); ++v) {
		if (_variables_is_fixed[v]) {
			if (_variables_type[v] == LinearSolver::FLOAT_VARIABLE)
				problem << " F_" << v ;
			else if (_variables_type[v] == LinearSolver::UINT_VARIABLE)
				problem << " UI_" << v ;
			else if (_variables_type[v] == LinearSolver::BIN_VARIABLE)
				problem << " B_" << v;
			else
				UNWANTED_BEHAVIOR("Unknown variable type, this is not supposed to happens, please investigate.");
			problem << " = " << _variables_value[v] << std::endl;
		}
	}

	problem << "Bounds " << std::endl;
	for (size_t v = 0; v < nb_variables(); ++v) if (_variables_bounds_index[v] != NOT_AN_ID) {
		if (_variables_bounds[_variables_bounds_index[v]][0] == LinearSolver::infinityminus)
			problem <<" -Inf <= ";
		else
			problem  << _variables_bounds[_variables_bounds_index[v]][0] << " <= ";
		if (_variables_type[v] == LinearSolver::FLOAT_VARIABLE)
			problem << " F_" << v;
		else if (_variables_type[v] == LinearSolver::UINT_VARIABLE)
			problem << " UI_" << v;
		else
			UNWANTED_BEHAVIOR("Unknown variable type, this is not supposed to happens, please investigate.");
		if (_variables_bounds[_variables_bounds_index[v]][1] == LinearSolver::infinityplus)
			problem << " <=  Inf " << std::endl;
		else
			problem << " <= " << _variables_bounds[_variables_bounds_index[v]][1] << std::endl;
	}
	problem << "General " << std::endl; // par default general est positif
	for (size_t v = 0; v < nb_variables(); ++v) {
		if (_variables_type[v] == LinearSolver::UINT_VARIABLE)
			problem << " UI_" << v;
	}
	problem << std::endl;


	bool has_binary = false;
	for (size_t v = 0; v < nb_variables(); ++v)
		has_binary = has_binary || (_variables_type[v] == LinearSolver::BIN_VARIABLE);
	if (has_binary) {
		problem << "Binary " << std::endl;
		for (size_t v = 0; v < nb_variables(); ++v) {
			if (_variables_type[v] == LinearSolver::BIN_VARIABLE)
				problem << " B_" << v;
		}
		problem << std::endl;
	}


	problem << "End " << std::endl;


	return problem.str();

}
bool LinearSolver::solve() {
	if (_is_in_constraint) UNWANTED_BEHAVIOR("End the constraint before solving");

	std::ofstream lp_file(_in_file);
	lp_file << print_problem();
	lp_file.close();
		
	//_out_file = "C:/NICO/prog/hexdom-with-geogram/test/tmp/solved.sol"; _is_solved=true;if(0)
	if (!run_solver()) {
		TALK("Problem while solving");
		return false;
	}
		
	if (!read_out_file()) {
		TALK("Problem while reading result, problem impossible ?");
		return false;
	}
	return true;
}


const std::string LinearSolver::print_solution() {
	if (!_is_solved) UNWANTED_BEHAVIOR("The problem was not solved, or modified after it was solved, can't print solution");
	std::ostringstream  solution;
	solution << "====== Printing solution ====== " << std::endl;;
	solution << print_timings();

	for (size_t v = 0; v < nb_variables(); ++v) {
		if (_variables_type[v] == LinearSolver::FLOAT_VARIABLE)
			solution << "F_" << v ;
		else if (_variables_type[v] == LinearSolver::UINT_VARIABLE)
			solution << "UI_" << v ;
		else if (_variables_type[v] == LinearSolver::BIN_VARIABLE)
			solution << "B_" << v;
		solution << "\t = " << _variables_value[v] << std::endl;
	}
	return solution.str();
}

const std::string LinearSolver::print_timings() {
	if (!_is_solved) UNWANTED_BEHAVIOR("The problem was not solved, or modified after it was solved, can't print solution");
	std::ostringstream  timing;
	timing << "Solution was found in " << _run_time << " sec" << std::endl;
	return timing.str();
}


double LinearSolver::get_variable(const size_t variable_id) {
	if (!_is_solved) UNWANTED_BEHAVIOR("The problem was not solved, or modified after it was solved, can't give a value");
	if (variable_id > nb_variables())
		UNWANTED_BEHAVIOR("out of bound of variable  : " << variable_id << " > " << nb_variables() - 1);
	return _variables_value[variable_id];
}

bool LinearSolver::run_solver() {
	bool went_well = true;
	const clock_t begin = clock();


	remove(_out_file.c_str());

	switch (_runned_solver) {
	case 1: // LP_SOLVE
	{
		std::string const cmd = lp_solve_binary_path + " -presolve -timeout " + std::to_string(_solver_time_out) + " -bfp bfp_LUSOL  -rxli xli_CPLEX " + _in_file + " > " + _out_file;
		TALK(cmd.c_str());
		if (std::system(cmd.c_str()) != 0) {
			TALK("Pb with LP_Solve");
			went_well = false;
		}
		TALK("solved.");
		break;
	}
	case 2: //CPLEX
	{

		std::string const cmd = cplex_binary_path + " -c \" read " + _in_file + " \" \"set threads 1 \" \" set timelimit " + std::to_string(_solver_time_out) +" \" \"opt \" \" set logfile " + _out_file + " \" \"display solution variables -\" \"quit\" > cplex.log ";
		TALK(cmd.c_str());
		//removing previous result file
		if (std::system(cmd.c_str()) != 0) {
			TALK("Pb with CPLEX");
			went_well = false;
		}
		// removing cplex logfile
		std::filesystem::remove("cplex.log");
		TALK("solved.");
		break;
	}
	case 3: //GPLK - GLPSOL
	{
		std::string const cmd = glpk_binary_path + " --tmlim " +std::to_string(_solver_time_out) +  "  --lp " + _in_file + " -o " + _out_file + " > glpk.log";
		TALK(cmd.c_str());
		if (std::system(cmd.c_str()) != 0) {
			TALK("Pb with glpsol");
			went_well = false;
		}
		std::filesystem::remove("glpk.log");

		TALK("solved.");
		break;
	}
	case 4: // GUROBI
	{
		std::string const cmd = gurobi_cl_binary_path + +" Threads=1 TimeLimit=" + std::to_string(_solver_time_out) + " ResultFile=" + _out_file + " " + _in_file + " > gurobi.log";
		TALK(cmd.c_str());
		//removing previous result file
		if (std::system(cmd.c_str()) != 0) {
			TALK("Pb with GUROBI");
			went_well = false;
		}
		// removing cplex logfile
		std::filesystem::remove("gurobi.log");
		TALK("solved.");

		break;
	}
	default:
		UNWANTED_BEHAVIOR("Value " << _runned_solver << " for solver type is unknown, stopping here. It wasn't suppose to got all the way here, please investigate.");
	}


	const clock_t end = clock();
	_run_time = (double)((int64_t)end - (int64_t)begin) / CLOCKS_PER_SEC;
	_is_solved = true;
	return went_well;
}
inline static bool string_start(std::string string, std::string start_of_string) {
	if (string.size() <= start_of_string.size()) return false;
	return (std::string(string.begin(), string.begin() + (long int)start_of_string.size()) == start_of_string);
}
bool LinearSolver::read_out_file() {
	TALK("Reading result...");
	std::ifstream solver_output(_out_file);
	std::string a;
	size_t compteur = 0;
	while (solver_output >> a && compteur < nb_variables()) {
		if (a.size() > 2)
		if (string_start(a, "F_") || (a.size() > 3 && string_start(a, "UI_")) || string_start(a, "B_")) {
			size_t id = 0;
			if (string_start(a, "F_")) {
				id = (size_t)stoi(std::string(a.begin() + 2, a.end()));
				//assert(_variables_type[id] == FLOAT_VARIABLE);
			}
			else if (string_start(a, "B_")) {
				id = (size_t)stoi(std::string(a.begin() + 2, a.end()));
				//assert(_variables_type[id] == BIN_VARIABLE);
				if (_runned_solver == 3) solver_output >> a; //necessary ?
			}
			else if (string_start(a, "UI_")) {
				id = (size_t)stoi(std::string(a.begin() + 3, a.end()));
				//assert(_variables_type[id] == UINT_VARIABLE);
				if (_runned_solver == 3) solver_output >> a;
			}
			while (compteur < id) {
				_variables_value[(compteur++)] = 0;
			}
			double tmp; solver_output >> tmp;
			_variables_value[(compteur++)] = tmp;//= index_t(stoi(a));
		}
	}

	solver_output.close();
	return true;

}



