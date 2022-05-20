#include <ultimaille/all.h>
#include "elliptic_smoothing.h"



bool smooth_tet_mesh(UM::Tetrahedra& m, std::vector<bool>& locks, const smoother_options& options);
bool smooth_tet_mesh(UM::Tetrahedra& m, const UM::Tetrahedra& ref, std::vector<bool>& locks, const smoother_options& options);

bool smooth_hex_mesh(UM::Hexahedra& m, const std::vector<bool>& locks);