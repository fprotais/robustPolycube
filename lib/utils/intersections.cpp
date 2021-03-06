#include "intersections.h"

using namespace UM; 

#define FOR(i, n) for(int i = 0; i < n; i++)

namespace intersections {
	bool point_is_in_tet(const vec3& A, const vec3& B, const vec3& C, const vec3& D, const vec3& P, std::array<double, 4>& l, const double dillatation) {
		constexpr int find[4][3] = { { 1, 2, 3 },{ 3, 2, 0 },{ 0, 1, 3 },{ 0, 2, 1 } };
		vec3 center = 0.25 * (A + B + C + D);
		double dillat = dillatation * (A - B).norm();
		center = center * dillat;
		std::array<vec3, 4> p = { { (1 + dillat) * A - center, (1 + dillat) * B - center, (1 + dillat) * C - center, (1 + dillat) * D - center } };

		double total_volume = (p[3] - p[0]) * cross(p[1] - p[0], p[2] - p[0]);

		FOR(f, 4) {
			const vec3 vA = p[find[f][0]];
			const vec3 vB = p[find[f][1]];
			const vec3 vC = p[find[f][2]];

			l[f] = (vC - P) * cross(vA - P, vB - P) / total_volume;
		}
		return (l[0] >= 0 && l[1] >= 0 && l[2] >= 0 && l[3] >= 0);

	}

	inline double det(const vec2& v1, const  vec2& v2) {
		return v1.x * v2.y - v1.y * v2.x;
	}
	bool point_is_in_triangle(const vec2& A, const vec2& B, const vec2& C, const vec2& P, std::array<double, 3>& l) {
		double area = det(B - A, C - A);
		l[0] = det(C - B, P - B) / area;
		l[1] = det(A - C, P - C) / area;
		l[2] = det(B - A, P - A) / area;

		return l[0]>= 0 && l[1]>=0 && l[2]>=0;
	}
	// code from : https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
	bool ray_trace_triangle(const vec3& A, const vec3& B, const vec3& C, const vec3& s1, const vec3& s2, double& l1, double& l2, double& l3) {
		vec3 E1 = B - A;
		vec3 E2 = C - A;
		vec3 N = cross(E1, E2);
		double det = -(s2 - s1) * N;
		double invdet = 1.0 / det;
		vec3 AO = s1 - A;
		vec3 DAO = cross(AO, s2 - s1);
		double u = E2 * DAO * invdet;
		double v = -E1 * DAO * invdet;
		double t = AO * N * invdet;
		l1 = (1 - u - v);
		l2 = u;
		l3 = v;
		return (det >= 1e-10 && t <= 1.0 && t >= 0.0 && u >= 0.0 && v >= 0.0 && (u + v) <= 1.0);
	}
	bool inter_sec_triangle(const vec3& A, const vec3& B, const vec3& C, const vec3& s1, const vec3& s2, double& l1, double& l2, double& l3) {
		if (ray_trace_triangle(A, B, C, s1, s2, l1, l2, l3))
			return true;
		else if (ray_trace_triangle(A, B, C, s2, s1, l1, l2, l3))
			return true;
		else
			return false;

	}
	bool inter_sec_quad(const vec3& A, const vec3& B, const vec3& C, const vec3& D, const vec3& s1, const vec3& s2, vec3& inter) {
		double l1, l2, l3;
		if (inter_sec_triangle(A, B, C, s1, s2, l1, l2, l3)) {
			inter = l1 * A + l2 * B + l3 * C;
			return true;
		} 
		if (inter_sec_triangle(C, D, A, s1, s2, l1, l2, l3)) {
			inter = l1 * C + l2 * D + l3 * A;
			return true;
		}
		return false;
	}

	bool point_is_cube(const vec3& O, const vec3& dO, const vec3& p) {
		FOR(dim, 3) if (p[dim]< O[dim] || p[dim]>O[dim] + dO[dim]) return false;
		return true;
	}
	static constexpr int facet_vertex[6][4] = { {0,4,6,2}, {1,3,7,5}, {0,1,5,4}, {2,6,7,3}, {0,2,3,1}, {4,5,7,6} };
	// check in 48 tets splitting the hex;
	bool point_is_in_hex(const std::array<UM::vec3, 8>& hex, const vec3& P) {
		vec3 cell_center;
		FOR(i, 8) cell_center += hex[i];
		cell_center = cell_center / 8;
		std::array<double, 4> l;
		FOR(cf, 6) {
			vec3 facet_center;
			FOR(i, 4) cell_center += hex[facet_vertex[cf][i]];
			facet_center = facet_center / 4;
			FOR(cfv, 4) {
				vec3 v1 = hex[facet_vertex[cf][cfv]];
				vec3 v2 = hex[facet_vertex[cf][(cfv + 1) % 4]];
				vec3 edge_center = 0.5 * (v1 + v2);
				if (point_is_in_tet(v1, edge_center, cell_center, facet_center, P, l)) return true;
				if (point_is_in_tet(edge_center, v2, cell_center, facet_center, P, l)) return true;
			}
		}
		return false;
	}

	bool point_is_in_quad(const std::array<UM::vec2, 4>& quad, const UM::vec2& P) {
		std::array<double, 3> l;
		if (point_is_in_triangle(quad[0], quad[1], quad[2], P, l)) return true;
		if (point_is_in_triangle(quad[0], quad[2], quad[3], P, l)) return true;
		return false;
	}



	bool inter_sec_sec(const UM::vec2& A1, const UM::vec2& B1, const UM::vec2& A2, const UM::vec2& B2, UM::vec2& P){
		vec2 s1 = A1 - B1;
		vec2 s2 = A2 - B2;
		vec2 s12 = A1 - A2;
		double d = det(s1, s2);
		if (d == 0) return false;

		double s = det(s12, s1) / d;
		double t = det(s12, s2) / d;
		if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
			P = A1 + t * s1;
			return true;
		}
		return false;
	}
}


