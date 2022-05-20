#include "pillowing.h"
#include <utils/standard_hexsmoother.h>

#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;

void add_a_pillow(Hexahedra & m, CellAttribute<int>&blockid) {
	std::vector<int> new_point(m.nverts(), -1);
	std::vector<std::array<int, 8>> new_hexes;
	OppositeFacet vec(m);
	int nb_new_points = 0;
	FOR(c, m.ncells()) FOR(cf, 6) if (vec.adjacent[6 * c + cf] == -1) {
		FOR(cfv, 4) if (new_point[m.facet_vert(c, cf, cfv)] == -1) {
			new_point[m.facet_vert(c, cf, cfv)] = m.nverts() + nb_new_points++;
		}
		std::array<int, 8> hex;
		FOR(cfv, 4) hex[cfv] = m.facet_vert(c, cf, cfv);
		FOR(cfv, 4) hex[cfv + 4] = new_point[m.facet_vert(c, cf, cfv)];
		std::swap(hex[2], hex[3]);
		std::swap(hex[6], hex[7]);
		new_hexes.push_back(hex);
	}
	int off_v = m.points.create_points(nb_new_points);
	FOR(v, off_v) if (new_point[v] != -1) m.points[new_point[v]] = m.points[v];
	int off_h = m.create_cells(new_hexes.size());
	FOR(i, new_hexes.size()) FOR(hv, 8) m.vert(off_h + i, hv) = new_hexes[i][hv];
	FOR(i, new_hexes.size()) blockid[i + off_h] = -1;


	OppositeFacet newvec(m);
	std::vector<bool> locks(m.nverts() * 3, false);
	FOR(c, m.ncells()) FOR(cf, 6) if (newvec.adjacent[6 * c + cf] == -1) FOR(cfv, 4) FOR(d, 3) locks[3 * m.facet_vert(c, cf, cfv) + d] = true;

	smooth_hex_mesh(m, locks);

}