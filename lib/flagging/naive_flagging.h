#pragma once
#include <ultimaille/all.h>

void generate_naive_flagging(const UM::Tetrahedra& m, UM::CellFacetAttribute<int>& flag);

void generate_naive_flagging(const UM::Triangles& m, UM::FacetAttribute<int>& flag);

