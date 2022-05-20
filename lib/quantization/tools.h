#pragma once
#include <ultimaille/all.h>
#include "dataStructure.h"





void disp_chart(const UM::Tetrahedra& m, const rb_data_structure::Chart& c, const std::string& name = "chart");

void disp_charts(const UM::Tetrahedra& m, const rb_data_structure::Sorted_charts& charts, const std::string& name = "charts");

void disp_charts_dim(const UM::Tetrahedra& m, const rb_data_structure::Sorted_charts& charts, int dim, const std::string& name = "charts");

void change_points_from_polycuboid_to_mesh(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, std::vector<UM::vec3>& points);


void disp_block_decomposition(rb_data_structure::Block_decomposition& blocks, const std::string& name = "blocks");
