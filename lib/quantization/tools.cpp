#include "tools.h"
#include <utils/trace.h>
#include <utils/projection.h>
#include <utils/intersections.h>

#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;

void disp_chart(const UM::Tetrahedra& m, const rb_data_structure::Chart& c, const std::string& name) {
	UM::CellFacetAttribute<int> chart(m, 0);
	for (int i : c.tetFacets) chart[i] = 1;
	Trace::drop_cellfacet_scalar(m, chart, name, 0, true);
}
void disp_charts(const UM::Tetrahedra& m, const rb_data_structure::Sorted_charts& charts, const std::string& name) {
	UM::CellFacetAttribute<int> chart(m, -1);
	FOR(cf, m.ncells() * 4) chart[cf] = charts.cf2chart[cf];
	Trace::drop_cellfacet_scalar(m, chart, name, -1, true);
}

void disp_charts_dim(const UM::Tetrahedra& m, const rb_data_structure::Sorted_charts& charts, int dim, const std::string& name) {
	UM::CellFacetAttribute<int> chart(m, -1);
	FOR(cf, m.ncells() * 4) if (charts.cf2chart[cf] != -1) if (charts.charts[charts.cf2chart[cf]].dim == dim) chart[cf] = charts.cf2chart[cf];
	Trace::drop_cellfacet_scalar(m, chart, name + "_dim_" + std::to_string(dim) + "_", -1, true);
}

void change_points_from_polycuboid_to_mesh(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, std::vector<vec3>& points) {
	std::vector<BBox3> boxes(polycuboid.ncells());
	FOR(c, polycuboid.ncells()) FOR(cv, 4) boxes[c].add(polycuboid.points[polycuboid.vert(c, cv)]);
	HBoxes matcher(boxes);

	for (vec3& P : points){
		BBox3 box; box.add(P);
		std::vector<int> candidates;
		matcher.intersect(box, candidates);

		std::vector<vec3> points_in_m;
		for (int c : candidates) {
			std::array<double, 4> l;
			vec3 A = polycuboid.points[polycuboid.vert(c, 0)];
			vec3 B = polycuboid.points[polycuboid.vert(c, 1)];
			vec3 C = polycuboid.points[polycuboid.vert(c, 2)];
			vec3 D = polycuboid.points[polycuboid.vert(c, 3)];
			if (intersections::point_is_in_tet(A, B, C, D, P, l)) {
				vec3 res;
				FOR(i, 4) res += l[i] * m.points[m.vert(c, i)];
				points_in_m.push_back(res);
			}
		}
		if (points_in_m.empty()){
			// no point in m, we project to closest
			double eps = 0.1;
			while (candidates.empty()) {
				box.dilate(eps);
				matcher.intersect(box, candidates);
				eps *= 2;
			}
			candidates.clear();
			matcher.intersect(box, candidates);

			double closest_dist = 1E100;
			vec3 closest_P;
			for (int c : candidates) FOR(cf, 4) {
				std::array<vec3, 3> ABC;
				FOR(i, 3) ABC[i] = m.points[m.facet_vert(c, cf, i)];
				std::array<double, 3> l;
				vec3 np;
				double dist = point_triangle_squared_distance(P, ABC, np, l);
				if (dist < closest_dist) {
					closest_dist = dist;
					closest_P = np;
				}
			}
			P = closest_P;
		}
		else {
			vec3 avg;
			for (vec3 np : points_in_m) avg += np;
			avg /= (double)points_in_m.size();
			P = avg;
		}
	}
}

void disp_block_decomposition(rb_data_structure::Block_decomposition& blocks, const std::string& name) {
	CellAttribute<int> id(blocks.m);
	FOR(c, blocks.m.ncells()) id[c] = c;
	Trace::drop_cells_scalar(blocks.m, id, name + "_on_mesh");
	
	FOR(v, blocks.m.nverts()) std::swap(blocks.m.points[v], blocks.block_coord[v]);
	Trace::drop_cells_scalar(blocks.m, blocks.cell_volume, name + "_blocks");
	FOR(v, blocks.m.nverts()) std::swap(blocks.m.points[v], blocks.block_coord[v]);

	FOR(v, blocks.m.nverts()) std::swap(blocks.m.points[v], blocks.float_coord[v]);
	Trace::drop_cells_scalar(blocks.m, id, name + "_on_polycuboid");
	FOR(v, blocks.m.nverts()) std::swap(blocks.m.points[v], blocks.float_coord[v]);

	FOR(v, blocks.m.nverts()) std::swap(blocks.m.points[v], blocks.int_coord[v]);
	Trace::drop_cells_scalar(blocks.m, id, name + "_on_polycube");
	FOR(v, blocks.m.nverts()) std::swap(blocks.m.points[v], blocks.int_coord[v]);
}
