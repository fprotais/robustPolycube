#pragma once
#include <ultimaille/all.h>
#include <array>
#include "dataStructure.h"

void simplify_block_structure(rb_data_structure::Block_decomposition& blocks);

void quantize(const rb_data_structure::Sorted_charts& charts, rb_data_structure::Block_decomposition& blocks, double scaling);

// quantize calls first MILPless_quantize before trying MILP solving
void MILPless_quantize(const rb_data_structure::Sorted_charts& charts, rb_data_structure::Block_decomposition& blocks, double scaling);

