#include "quantization.h"
#include <algorithm>
#include <utils/trace.h>
#include <utils/LinearSolver.h>

#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;
using namespace rb_data_structure;

void simplify_block_structure(rb_data_structure::Block_decomposition & blocks) {
	Trace::alert("Pouet, function simplify_block_structure need implementation.");
}

// fixing 3 variables quite far : (invariant by translation, one in each direction)
static void fixe_0(LinearSolver& solver, const Sorted_charts& charts, double far_shift) {
	solver.fix_variable(charts.ranked(0, 0), far_shift);
	solver.fix_variable(charts.ranked(1, 0), far_shift);
	solver.fix_variable(charts.ranked(2, 0), far_shift);
}

static void initialise_energy_on_blocks(LinearSolver& solver, const Sorted_charts& charts, const Block_decomposition& blocks, double scale) {
	FOR(c, blocks.m.ncells()) FOR(d, 3) {
		int alpha_id = solver.add_variable(blocks.cell_volume[c], LinearSolver::FLOAT_VARIABLE);
		int c1 = blocks.issuing_chart_id[6 * c + 2 * d + 0];
		int c2 = blocks.issuing_chart_id[6 * c + 2 * d + 1];
		if (charts[c1].value < charts[c2].value) std::swap(c1, c2);
		//if (is_shadow_chart[c1] || is_shadow_chart[c2]) continue;
		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(alpha_id, 1);
		solver.add_constraint_coefficient(c1, -1);
		solver.add_constraint_coefficient(c2, 1);
		solver.add_constraint_righthandside(-scale * (charts[c1].value - charts[c2].value));
		solver.end_constraint();
		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(alpha_id, 1);
		solver.add_constraint_coefficient(c1, 1);
		solver.add_constraint_coefficient(c2, -1);
		solver.add_constraint_righthandside(scale * (charts[c1].value - charts[c2].value));
		solver.end_constraint();
	}
}

static std::vector<int> compute_vert_connectivity_graph(const Block_decomposition& blocks) {
	std::vector<int> graph(blocks.m.nverts() * 6, -1);

	FOR(c, blocks.m.ncells()) FOR(cf, 6) FOR(cfv, 4) {
		int v1 = blocks.m.facet_vert(c, cf, cfv);
		int v2 = blocks.m.facet_vert(c, cf, (cfv + 1) % 4);
		FOR(dim, 3) if (blocks.block_coord[v1][dim] > blocks.block_coord[v2][dim] + 0.5) {
			graph[6 * v2 + dim * 2 + 0] = v1;
			graph[6 * v1 + dim * 2 + 1] = v2;
		}
	}
	return graph;
}

static void initialise_straight_path_constraint(LinearSolver& solver, const Sorted_charts& charts, const Block_decomposition& blocks, const std::vector<int>& graph) {
	std::vector<std::array<int, 3>> vert_charts(blocks.m.nverts(), { -1,-1,-1 });
	
	FOR(c, blocks.m.ncells())  FOR(cf, 6) if (blocks.bnd_chart_id[6 * c + cf] != -1) {
		int ch = blocks.bnd_chart_id[6 * c + cf];
		int d = charts[ch].dim;
		FOR(cfv, 4) vert_charts[blocks.m.facet_vert(c, cf, cfv)][d] = ch;
	}

	FOR(v, blocks.m.nverts()) FOR(d, 3) if (vert_charts[v][d] != -1) {
		int direction = 2 * d;
		int nv = graph[6 * v + direction];
		if (nv == -1) continue; // only making upward path to avoir repetition in constraints
		while (vert_charts[nv][d] == -1) {
			nv = graph[6 * nv + direction];
			if (nv == -1) {
				Trace::alert("Should not happen?");
				abort();
			}
		}
		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(vert_charts[nv][d], 1);
		solver.add_constraint_coefficient(vert_charts[v][d], -1);
		solver.add_constraint_righthandside(1);
		solver.end_constraint();
	}
}

static void initialise_global_order_constraint(LinearSolver& solver, const Sorted_charts& charts, const Block_decomposition& blocks) {
	FOR(dim, 3) FOR(c, charts.nb_charts(dim) - 1) {
		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(charts.ranked(dim, c + 1), 1);
		solver.add_constraint_coefficient(charts.ranked(dim, c), -1);
		solver.end_constraint();
	}
}
static void initialise_local_order_constraint(LinearSolver& solver, const Sorted_charts& charts, const Block_decomposition& blocks) {
	FOR(c, blocks.m.ncells()) FOR(d, 3) {
		int alpha_id = solver.add_variable(blocks.cell_volume[c], LinearSolver::FLOAT_VARIABLE);
		int c1 = blocks.issuing_chart_id[6 * c + 2 * d + 0];
		int c2 = blocks.issuing_chart_id[6 * c + 2 * d + 1];
		int dim = charts[c1].dim;
		if (charts[c1].value < charts[c2].value) std::swap(c1, c2); // we could compare c1 and c2, they are sorted...
		solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
		solver.add_constraint_coefficient(c1, 1);
		solver.add_constraint_coefficient(c2, -1);
		solver.end_constraint();
		
	}
}

static bool check_path_constraint_in_R(const Sorted_charts& charts, const Block_decomposition& blocks, const std::vector<int>& chart_value, std::array<int, 6>& new_constraint) {
	std::vector<std::vector<int>> chart_verts(charts.size());
	FOR(c, blocks.m.ncells())  FOR(cf, 6) if (blocks.bnd_chart_id[6*c+cf] != -1) {
		FOR(cfv, 4) chart_verts[blocks.bnd_chart_id[6 * c + cf]].push_back(blocks.m.facet_vert(c, cf, cfv));
	}
	FOR(c, chart_verts.size()) {
		std::sort(chart_verts[c].begin(), chart_verts[c].end());
		auto last = std::unique(chart_verts[c].begin(), chart_verts[c].end());
		chart_verts[c].erase(last, chart_verts[c].end());
	}
	FOR(dim, 3) {
		FOR(i1, charts.nb_charts(dim) - 1) for (int i2 = charts.nb_charts(dim) - 1; i2 > i1; i2--){
			int c1 = charts.ranked(dim, i1);
			int c2 = charts.ranked(dim, i2);
			if (chart_value[c1] != chart_value[c2])  continue;
			//if (is_shadow_chart[c1] || is_shadow_chart[c2]) continue;
			for (int v1 : chart_verts[c1]) for (int v2 : chart_verts[c2]) {
				int C1_[3], C2_[3];
				FOR(d2, 3) {
					C1_[d2] = charts.ranked(d2, (int)blocks.block_coord[v1][d2]);
					C2_[d2] = charts.ranked(d2, (int)blocks.block_coord[v2][d2]);
				}

				if (chart_value[C1_[0]] == chart_value[C2_[0]]
					&& (chart_value[C1_[1]] == chart_value[C2_[1]])
					&& (chart_value[C1_[2]] == chart_value[C2_[2]])) {
					new_constraint = { C1_[0], C2_[0], C1_[1], C2_[1], C1_[2], C2_[2] };
					FOR(d, 3) if (new_constraint[2 * d] < new_constraint[2 * d + 1]) std::swap(new_constraint[2 * d], new_constraint[2 * d + 1]);
					std::cerr << "Two charts were merged, adding a constraint : " << c1 << "," << c2 << " (dim = " << dim << ")" << std::endl;;
					return false;
				}
			}

		}
	}
	return true;
}

static bool check_path_constraint_in_blocks(const Sorted_charts& charts, const Block_decomposition& blocks, const std::vector<int>& chart_values, const std::vector<int>& graph, std::array<int, 6>& new_constraint) {
	Trace::alert("Function check_path_constraint_in_blocks is not implemented");
	return true;
}

static void add_constraint(LinearSolver& solver, const std::array<int, 6>& new_constraint) {
	solver.begin_constraint(LinearSolver::INEQUALITY_SUP);
	FOR(d, 3) {
		solver.add_constraint_coefficient(new_constraint[2*d], 1);
		solver.add_constraint_coefficient(new_constraint[2*d + 1], -1);
		
	}
	solver.add_constraint_righthandside(1);
	solver.end_constraint();
}

void quantize(const rb_data_structure::Sorted_charts & charts, rb_data_structure::Block_decomposition& blocks, double scale) {
	FOR(d, 3) if (charts.nb_charts(d) == 0) {
		Trace::alert("There is no chart in dimension " + std::to_string(d) + ". This is code breaking, we abort.");
		abort();
	}

	MILPless_quantize(charts, blocks, scale);

	if (!LinearSolver_is_runnable) {
		Trace::warning("no available MILP solver. Only MILPless_quantize was run. See README.md on how to set-up MILP solvers (If I had time to fill it, or just contact author).");
		return;
	}

	LinearSolver solver(LinearSolver::MINIMIZE, 40);
	solver.add_variables(std::vector<double>(charts.size(), 0), std::vector<LinearSolver::variable_types>(charts.size(), LinearSolver::UINT_VARIABLE));
	
	const int far_shift = std::max(500, (int)charts.size());

	fixe_0(solver, charts, far_shift);

	std::vector<int> vert_connectivity_graph = compute_vert_connectivity_graph(blocks);

	initialise_energy_on_blocks(solver, charts, blocks, scale);



	const bool allow_overlaps = false;

	if (!allow_overlaps)
		initialise_global_order_constraint(solver, charts, blocks);
	else
		initialise_local_order_constraint(solver, charts, blocks);

	initialise_straight_path_constraint(solver, charts, blocks, vert_connectivity_graph);

	std::vector<int> chart_values(charts.size(), 0);
	bool is_valid = false;
	int nb_of_iter = 0;
	while (!is_valid && nb_of_iter++ < 100) {
		std::cerr << "============= ITER : " << nb_of_iter << " =============" << std::endl;
		solver.solve();
		std::cerr << solver.print_timings() << std::endl;

		FOR(c, charts.size()) chart_values[c] = int(solver.get_variable(c));
		std::array<int, 6> new_constraint;
		if (!allow_overlaps)
			is_valid = check_path_constraint_in_R(charts, blocks, chart_values, new_constraint);
		else
			is_valid = check_path_constraint_in_blocks(charts, blocks, chart_values, vert_connectivity_graph, new_constraint);
		if (!is_valid) add_constraint(solver, new_constraint);
	}
	if (nb_of_iter >= 100) {
		Trace::alert("We iterated to much, no valid solution was found, stopping, a solution was still found earlier.");
	}
	if (is_valid) {
		FOR(c, blocks.m.ncells()) FOR(cf, 6) {
			int ch = blocks.issuing_chart_id[6 * c + cf];
			int dim = charts[ch].dim;
			FOR(cfv, 4) blocks.int_coord[blocks.m.facet_vert(c, cf, cfv)][dim] = chart_values[ch] - far_shift;
		}
	}

}



void MILPless_quantize(const rb_data_structure::Sorted_charts& charts, rb_data_structure::Block_decomposition& blocks, double scaling) {
	FOR(d, 3) if (charts.nb_charts(d) == 0) {
		Trace::alert("There is no chart in dimension " + std::to_string(d) + ". This is code breaking, we abort.");
		abort();
	}

	std::vector<int> init_chart_values(charts.size(), 0);
	FOR(c, charts.size()) init_chart_values[c] = (int)std::round(scaling * charts[c].value);

	std::vector<int> vert_connectivity_graph = compute_vert_connectivity_graph(blocks);

	std::vector<int> chart_values = init_chart_values;
	const bool allow_overlaps = false;
	std::cerr << "MILPless quantization:" << std::endl;
	bool is_valid = false;
	int nb_of_iter = 0;
	while (!is_valid && nb_of_iter++ < 100) {
		std::array<int, 6> new_constraint;
		if (!allow_overlaps)
			is_valid = check_path_constraint_in_R(charts, blocks, chart_values, new_constraint);
		else
			is_valid = check_path_constraint_in_blocks(charts, blocks, chart_values, vert_connectivity_graph, new_constraint);
		
		if (is_valid) break;
		
		int max_d = 0;
		double min_val = 1E100;
		FOR(d, 3) {
			int c1 = new_constraint[d * 2];
			int c2 = new_constraint[d * 2 + 1];
			double deltag = charts[c1].value - charts[c2].value;
			double deltaint = chart_values[c1] - chart_values[c2];
			double val = std::abs(1 - (deltaint + 1) / deltag);
			if (min_val > val) {
				max_d = d;
				min_val = val;
			}
		}
		for (int c_ = new_constraint[2 * max_d]; c_ < charts.start_of_dim[max_d] + charts.nb_charts(max_d); c_++) chart_values[c_]++;
	}
	if (nb_of_iter >= 100) {
		FOR(c, charts.size()) chart_values[c] = init_chart_values[c] + charts.dim_rank(c);
	}
	FOR(c, blocks.m.ncells()) FOR(cf, 6) {
		int ch = blocks.issuing_chart_id[6 * c + cf];
		int dim = charts[ch].dim;
		FOR(cfv, 4) blocks.int_coord[blocks.m.facet_vert(c, cf, cfv)][dim] = chart_values[ch];
	}
}