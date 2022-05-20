#include "naive_flagging.h"

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

const std::array<UM::vec3, 6> flag2normal = { UM::vec3(1,0,0), UM::vec3(-1,0,0), UM::vec3(0,1,0), UM::vec3(0,-1,0), UM::vec3(0,0,1), UM::vec3(0,0,-1) };


void generate_naive_flagging(const UM::Tetrahedra& m, UM::CellFacetAttribute<int>& flag) {
	OppositeFacet vec(m);
	FOR(cf, m.ncells() * 4) if (vec.adjacent[cf] == -1) {
		vec3 n = m.util.facet_normal(cf / 4, cf % 4);
		flag[cf] = 0;
		double best = flag2normal[0] * n;
		FOR(i, 5) {
			double d = flag2normal[i + 1] * n;
			if (d > best) {
				best = d;
				flag[cf] = i + 1;
			}
		}
	}
	else flag[cf] = -1;
}

void generate_naive_flagging(const UM::Triangles& m, UM::FacetAttribute<int>& flag) {
	FOR(f, m.nfacets()) {
		vec3 n = m.util.normal(f);
		flag[f] = 0;
		double best = flag2normal[0] * n;
		FOR(i, 5) {
			double d = flag2normal[i + 1] * n;
			if (d > best) {
				best = d;
				flag[f] = i + 1;
			}
		}
	}
}