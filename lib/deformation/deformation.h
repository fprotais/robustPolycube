#pragma once
#include <ultimaille/all.h>

// if constrained_chart_values is true, we lock chart's points to the position they are in U
void cube_cover(const UM::Tetrahedra& m, const UM::CellFacetAttribute<int>& flag, UM::PointAttribute<UM::vec3>& U, bool constrained_charts_values = false);

void correct_param(const UM::Tetrahedra& m, const UM::CellFacetAttribute<int>& flag, UM::PointAttribute<UM::vec3>& U, bool constrained_charts_values = false);
