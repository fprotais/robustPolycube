#pragma once
#include <ultimaille/all.h>
#include <array>
#include "dataStructure.h"

bool polycuboid_is_valid(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, UM::CellFacetAttribute<int>& cfflag);

bool cleanflags(const UM::Tetrahedra& m, const UM::Tetrahedra& polycuboid, UM::CellFacetAttribute<int>& cfflag);
