#include "tools.h"

#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;

double avgedgesize(const UM::Hexahedra& m){
 double avg = 0;
 FOR(c, m.ncells()) FOR(cf, 6) FOR(cfv, 4) {
	 avg += (m.points[m.facet_vert(c, cf, cfv)] - m.points[m.facet_vert(c, cf, (cfv + 1) % 4)]).norm();
 }
 avg /= 24 * m.ncells();
 return avg;
}