#pragma once
#include <ultimaille/all.h>
#include <array>
#include "dataStructure.h"

void inverse(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, const rb_data_structure::Sorted_charts& charts, const rb_data_structure::Block_decomposition& blocks, rb_data_structure::Final_mesh& Final_mesh);