#include "tools.h"
#include <utils/trace.h>
#include <fstream>
#include <array>

#include <utils/intersections.h>

#include "tetgen.h"


using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)



void disp_polycube(UM::Tetrahedra& m, PointAttribute<vec3>& U, UM::CellFacetAttribute<int>& cfflags, const std::string& fn) {
	FOR(v, m.nverts()) std::swap(m.points[v], U[v]);
	Trace::drop_volume(m, fn + "_polycuboid", {});
	Trace::drop_cellfacet_scalar(m, cfflags, "_wflagging", -1, true);
	FOR(v, m.nverts()) std::swap(m.points[v], U[v]);
	Trace::drop_volume(m, fn + "_param", { {{"U", U.ptr}},{},{},{}});
}




void transfer_to_surf(const UM::Tetrahedra& vol, const UM::CellFacetAttribute<int>& cflags, UM::Triangles& m, UM::FacetAttribute<int>& flags) {
	{
		std::vector<bool> del(m.nfacets(), true);
		m.delete_facets(del);
		m.delete_isolated_vertices();
	}
	OppositeFacet vec(vol);
	m.points.data->assign(vol.points.begin(), vol.points.end());
	FOR(cf, 4 * vol.ncells()) if (cflags[cf] != -1) {
		int f = m.create_facets(1);
		flags[f] = cflags[cf];
		FOR(fv, 3) m.vert(f, fv) = vol.facet_vert(cf / 4, cf % 4, fv);
	}
	m.delete_isolated_vertices();
}




static void add_bnd_quads(Triangles& m, FacetAttribute<int>& flag, vec3& point_out) {
	//if (m.nfacets() > 0) {
	//	vec3 n = cross(m.points[m.vert(0, 1)] - m.points[m.vert(0, 0)], m.points[m.vert(0, 2)] - m.points[m.vert(0, 0)]);
	//	vec3 bary = m.points[m.vert(0, 0)] + m.points[m.vert(0, 1)] + m.points[m.vert(0, 2)];
	//	bary /= 3;
	//	point_in = bary - 0.01 * n;
	//}
	//else {
	//	point_in = { -1,-1,-1 };
	//}


	vec3 bboxmin, bboxmax;
	bboxmin = m.points[0];
	bboxmax = m.points[0];
	FOR(v, m.nverts()) FOR(d, 3) {
		bboxmax[d] = std::max(bboxmax[d], m.points[v][d]);
		bboxmin[d] = std::min(bboxmin[d], m.points[v][d]);
	}
	vec3 diag = (bboxmax - bboxmin);
	bboxmax = bboxmax + 4*diag;
	bboxmin = bboxmin - 4*diag;
	point_out = bboxmin + 0.1 * diag;

	int off_v = m.points.create_points(14);
	m.points[off_v + 0] = bboxmin;
	m.points[off_v + 1] = { bboxmax[0], bboxmin[1], bboxmin[2] };
	m.points[off_v + 2] = { bboxmin[0], bboxmax[1], bboxmin[2] };
	m.points[off_v + 3] = { bboxmax[0], bboxmax[1], bboxmin[2] };
	m.points[off_v + 4] = { bboxmin[0], bboxmin[1], bboxmax[2] };
	m.points[off_v + 5] = { bboxmax[0], bboxmin[1], bboxmax[2] };
	m.points[off_v + 6] = { bboxmin[0], bboxmax[1], bboxmax[2] };
	m.points[off_v + 7] = bboxmax;


	m.points[off_v + 8 + 3] = { bboxmin[0], 0.5 * (bboxmin[1] + bboxmax[1]), 0.5 * (bboxmin[2] + bboxmax[2]) };
	m.points[off_v + 8 + 2] = { bboxmax[0], 0.5 * (bboxmin[1] + bboxmax[1]), 0.5 * (bboxmin[2] + bboxmax[2]) };

	m.points[off_v + 8 + 0] = { 0.5 * (bboxmin[0] + bboxmax[0]), bboxmin[1], 0.5 * (bboxmin[2] + bboxmax[2]) };
	m.points[off_v + 8 + 4] = { 0.5 * (bboxmin[0] + bboxmax[0]), bboxmax[1], 0.5 * (bboxmin[2] + bboxmax[2]) };

	m.points[off_v + 8 + 5] = { 0.5 * (bboxmin[0] + bboxmax[0]), 0.5 * (bboxmin[1] + bboxmax[1]), bboxmin[2] };
	m.points[off_v + 8 + 1] = { 0.5 * (bboxmin[0] + bboxmax[0]), 0.5 * (bboxmin[1] + bboxmax[1]), bboxmax[2] };




	int off_f = m.create_facets(24);
	int quads[6][4] = {
		{0,1,5,4},
		{4,5,7,6},
		{1,3,7,5},
		{0,4,6,2},
		{2,6,7,3},
		{0,2,3,1},
	};
	int quadflag[6] = { 3,4,0,1,2,5 };
	FOR(i, 6) {
		FOR(j, 4) {
			m.vert(off_f + 4 * i + j, 0) = off_v + quads[i][j];
			m.vert(off_f + 4 * i + j, 1) = off_v + quads[i][(j + 1) % 4];
			m.vert(off_f + 4 * i + j, 2) = off_v + 8 + i;
			flag[off_f + 4 * i + j] = quadflag[i];
		}
	}
}

inline std::array<int, 3> ordered_array(const std::array<int, 3>& arr) {
	std::array<int, 3> res = arr;
	if (res[0] > res[1]) std::swap(res[0], res[1]);
	if (res[1] > res[2]) std::swap(res[1], res[2]);
	if (res[0] > res[1]) std::swap(res[0], res[1]);
	return res;
}

inline bool equal(const std::array<int, 3>& arr1, const std::array<int, 3>& arr2) {
	return (arr1[0] == arr2[0] && arr1[1] == arr2[1] && arr1[2] == arr2[2]);
}
static void put_in_box(UM::Triangles& m, UM::FacetAttribute<int>& flag, const vec3& point_in, UM::Tetrahedra& box, UM::CellFacetAttribute<int>& box_flag, UM::CellAttribute<int>& is_in, int verbose) {
	vec3 point_out;
	add_bnd_quads(m, flag, point_out);
	if (verbose > 1) Trace::drop_facet_scalar(m, flag, "boxsurf_wflag");

	tetgenio tetgen_out_;
	tetgenio tetgen_in_;
	tetgenbehavior tetgen_args_;
	tetgen_in_.initialize();
	tetgen_in_.firstnumber = 0;
	tetgen_in_.pointlist = nullptr;
	if (verbose) tetgen_args_.parse_commandline("Vpq1.2/18YAAnn");
	else tetgen_args_.parse_commandline("Qpq1.2/18YAAnn");


	tetgen_in_.numberofpoints = m.nverts();
	tetgen_in_.pointlist = new double[3 * tetgen_in_.numberofpoints];
	FOR(v, m.nverts()) FOR(d, 3) tetgen_in_.pointlist[3 * v + d] = m.points[v][d];

	tetgenio::polygon* polygons = new tetgenio::polygon[m.nfacets()];
	tetgen_in_.numberoffacets = int(m.nfacets());
	tetgen_in_.facetmarkerlist = new int[m.nfacets()];
	tetgen_in_.facetlist = new tetgenio::facet[tetgen_in_.numberoffacets];
	FOR(f, m.nfacets()) {
		tetgenio::facet& F = tetgen_in_.facetlist[f];
		tetgenio::init(&F);
		F.numberofpolygons = 1;
		F.polygonlist = &polygons[f];
		tetgenio::polygon& P = F.polygonlist[0];
		tetgenio::init(&P);
		P.numberofvertices = int(m.facet_size(f));
		P.vertexlist = new int[m.facet_size(f)];
		FOR(fv, m.facet_size(f)) P.vertexlist[fv] = m.vert(f, fv);
		F.numberofholes = 0;
		F.holelist = nullptr;
		tetgen_in_.facetmarkerlist[f] = flag[f] + 1;
	}
	tetgen_in_.numberofregions = 2;
	tetgen_in_.regionlist = new double[10];
	FOR(d, 3) tetgen_in_.regionlist[d] = point_in[d];
	tetgen_in_.regionlist[3] = 1;
	tetgen_in_.regionlist[4] = 1.;
	FOR(d, 3) tetgen_in_.regionlist[5 + d] = point_out[d];
	tetgen_in_.regionlist[8] = 0;
	tetgen_in_.regionlist[9] = 1.;

	bool there_was_an_error = false;
	int error_code = 0;
	tetrahedralize(
		&tetgen_args_, &tetgen_in_, &tetgen_out_
	);


	// Deallocate the datastructures used by tetgen,
	// and disconnect them from tetgen,
	// so that tetgen does not try to deallocate them.

	// Pointlist was allocated in local array
	tetgen_in_.numberofpoints = 0;
	delete[] tetgen_in_.pointlist;
	tetgen_in_.pointlist = nullptr;

	// Edges were shared with constraint mesh
	// (no need to deallocate)
	tetgen_in_.numberofedges = 0;
	tetgen_in_.edgelist = nullptr;

	// Facets structures were allocated in local
	// array, and vertices indices were shared
	// with constraint mesh
	delete[] tetgen_in_.facetlist;
	delete[] tetgen_in_.facetmarkerlist;
	tetgen_in_.facetlist = nullptr;
	tetgen_in_.facetmarkerlist = nullptr;
	tetgen_in_.numberoffacets = 0;
	delete[] polygons;

	tetgen_in_.numberofregions = 0;
	delete[] tetgen_in_.regionlist;
	tetgen_in_.regionlist = nullptr;

	box.points.create_points(tetgen_out_.numberofpoints);
	FOR(v, box.nverts()) FOR(d, 3) box.points[v][d] = tetgen_out_.pointlist[3 * v + d];

	box.create_cells(tetgen_out_.numberoftetrahedra);
	FOR(tet, box.ncells()) {
		FOR(tet_vert, 4) box.vert(tet, tet_vert) = tetgen_out_.tetrahedronlist[4 * tet + tet_vert];
		is_in[tet] = (int)tetgen_out_.tetrahedronattributelist[tet];
		FOR(cf, 4) box_flag[4 * tet + cf] = -1;
	}
	FOR(tri, tetgen_out_.numberoftrifaces) {
		int f = tetgen_out_.trifacemarkerlist[tri] - 1;
		std::array<int, 3> v1;
		FOR(i, 3) v1[i] = tetgen_out_.trifacelist[3 * tri + i];
		FOR(j, 2) {
			int c = tetgen_out_.face2tetlist[2 * tri + j];
			if (c == -1) continue;
			FOR(cf, 4) {
				std::array<int, 3> v2;
				FOR(cfv, 3) v2[cfv] = box.facet_vert(c, cf, cfv);
				if (equal(ordered_array(v1), ordered_array(v2))) {
					box_flag[4 * c + cf] = f;
				}
			}
		}
	}
	OppositeFacet vec(box);
	FOR(c, box.ncells()) if (!is_in[c]) FOR(cf, 4) if (vec.adjacent[4 * c + cf] != -1) {
		box_flag[4 * c + cf] = -1;
	}
	if (verbose > 1) Trace::drop_cells_scalar(box, is_in, "box_in", 0);
	if (verbose > 1) Trace::drop_cells_scalar(box, is_in, "box_out", 1);
}



void put_model_in_box(const UM::Triangles& input, const UM::FacetAttribute<int>& flag, UM::Tetrahedra& m, UM::CellFacetAttribute<int>& cfflags, UM::CellAttribute<int>& inmesh) {
	std::vector<bool> del(m.ncells(), true);
	m.delete_cells(del);
	m.delete_isolated_vertices();

	vec3 n = cross(input.points[input.vert(0, 1)] - input.points[input.vert(0, 0)], input.points[input.vert(0, 2)] - input.points[input.vert(0, 0)]);
	vec3 point_near_bnd = input.util.bary_verts(0) - 0.01 * n;


	// manual copy of input to allow modification in put_in_box
	Triangles surf;
	surf.points.data->assign(input.points.begin(), input.points.end());
	surf.facets = input.facets;
	FacetAttribute<int> surf_flags(surf);
	FOR(f, surf.nfacets()) surf_flags[f] = flag[f];
	std::cerr << "Calling Tetgen..." << std::endl;
	put_in_box(surf, surf_flags, point_near_bnd, m, cfflags, inmesh, 0);
}


void remove_outerbox(UM::Tetrahedra& m, UM::CellFacetAttribute<int>& flags, UM::CellAttribute<int>& inmesh) {
	std::vector<bool> del(m.ncells());
	FOR(c, m.ncells()) del[c] = !inmesh[c];
	m.delete_cells(del);
	m.delete_isolated_vertices();
}



double avg_edge_size(const Triangles& m) {
	double size = 0;
	FOR(f, m.nfacets()) FOR(fc, 3) {
		size += (m.points[m.vert(f, fc)] - m.points[m.vert(f, (fc + 1) % 3)]).norm();
	}
	size /= m.nfacets() * 3;
	return size;
}

std::tuple<double, UM::vec3>  center_and_normalise_mesh(Triangles& m) {
	double scale = avg_edge_size(m);
	vec3 center;
	FOR(v, m.nverts()) center += m.points[v];
	center /= m.nverts();
	FOR(v, m.nverts()) m.points[v] = (m.points[v] - center) / scale;
	return {scale, center};
}

