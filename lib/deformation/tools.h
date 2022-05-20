#pragma once
#include <ultimaille/all.h>
#include <array>


void disp_polycube(UM::Tetrahedra& m, UM::PointAttribute<UM::vec3>& U, UM::CellFacetAttribute<int>& cfflags, const std::string& fn);

void transfer_to_surf(const UM::Tetrahedra& vol, const UM::CellFacetAttribute<int>& cflags, UM::Triangles& m, UM::FacetAttribute<int>& flags);

void put_model_in_box(const UM::Triangles& input, const UM::FacetAttribute<int>& flag, UM::Tetrahedra& m, UM::CellFacetAttribute<int>& cfflags, UM::CellAttribute<int>& inmesh);

void remove_outerbox(UM::Tetrahedra& model, UM::CellFacetAttribute<int>& flags, UM::CellAttribute<int>& inmesh);

double avg_edge_size(const UM::Triangles& m);

void center_and_normalise_mesh(UM::Triangles& m);