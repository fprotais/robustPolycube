#pragma once
#include <ultimaille/all.h>
#include <array>
#include "dataStructure.h"

void compute_charts(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, const UM::CellFacetAttribute<int>& cfflag, rb_data_structure::Sorted_charts& charts);

void chart_voxelisation(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, const rb_data_structure::Sorted_charts& charts, rb_data_structure::Block_decomposition& blocks);

