#ifndef __TRACE_H__
#define __TRACE_H__




#include <ultimaille/all.h>
#include <iostream>
#include <string>

#include <fstream>
#include <math.h>

#if __cplusplus >= 199711L  
#include <filesystem>
#endif

double ave_edge_size(const UM::Surface& m);

namespace Trace {

	extern bool drop_mesh_is_active;

	//    _      _ _   _      _ _         _   _          
	//   (_)_ _ (_) |_(_)__ _| (_)_____ _| |_(_)___ _ __ 
	//   | | ' \| |  _| / _` | | |_ / _` |  _| / _ \ '  | 
	//   |_|_||_|_|\__|_\__,_|_|_/__\__,_|\__|_\___/_||_|
	void initialize(const std::string& graphite_path);
	void conclude();

	//    _                _               _   ___ ___ 
	//   | |   ___  __ _  | |_ _____ __   /_\ | _ \_ _|
	//   | |__/ _ \/ _` | |  _/ -_) \ /  / _ \|  _/| | 
	//   |____\___/\__, |  \__\___/_\_\ /_/ \_\_| |___|
	//             |___/                              
	struct SwitchTextInScope {
		SwitchTextInScope(bool make_active);
		~SwitchTextInScope();
		static void force_activity(bool b); // dirty, doesn't stack change
		bool save_val;
	};


	struct Section {
		Section(const std::string& str);
		~Section();
		bool trace_was_active;
	};
	void step(std::string stepname, int alert_level = 0);
	void alert(std::string msg);

	void log_value(std::string const& str, double val, int alert_level = 0);
	void log_string(std::string const& str, std::string const& val, int alert_level = 0);
	void show_log();
	void append_py_log(std::string const& filename);


	//    ___                 __  __        _               _   ___ ___ 
	//   |   \ _ _ ___ _ __  |  \/  |___ __| |_  ___ ___   /_\ | _ \_ _|
	//   | |) | '_/ _ \ '_ \ | |\/| / -_|_-< ' \/ -_|_-<  / _ \|  _/| | 
	//   |___/|_| \___/ .__/ |_|  |_\___/__/_||_\___/__/ /_/ \_\_| |___|
	//                |_|                                               
	struct SwitchDropInScope {
		SwitchDropInScope(bool make_active);
		~SwitchDropInScope();
		bool save_val;
		static void force_activity(bool b); // dirty, doesn't stack change
	};


	void drop_volume(const UM::Volume& m,
		const std::string& name,
		const UM::VolumeAttributes& attr,
		std::string cell_attr = "",
		std::string point_attr = ""
	);

	void drop_surface(const UM::Surface& m,
		const std::string& name,
		const UM::SurfaceAttributes& attr,
		std::string facet_attr = "",
		std::string point_attr = "",
		std::string corner_attr = ""
	);

	void drop_polyline(const UM::PolyLine& m,
		const std::string& name,
		std::vector<UM::NamedContainer> vertex_attributes,
		std::vector<UM::NamedContainer> segment_attributes,
		bool show_edge_attr = false);


	struct ArrowStyle {
		ArrowStyle(int p_resolution, double p_diameter = 0) { resolution = p_resolution; diameter = p_diameter; }
		int resolution = 3;
		double diameter = 0;
	};

	void drop_facet_vec3(const UM::Surface& m, UM::FacetAttribute<UM::vec3>& facet_attr, std::string name = "vec3", ArrowStyle as = ArrowStyle(0, 0), bool proportional_to_ave_edge_size = true, bool normalize = false);


	void drop_surface_points_vec3(const UM::Surface& m, UM::PointAttribute<UM::vec3>& attr, std::string name = "vec3", ArrowStyle as = ArrowStyle(0, 0), bool proportional_to_ave_edge_size = true, bool normalize = false);

	void drop_texture(const UM::Surface& m, const std::string& name, UM::CornerAttribute<UM::vec2>& u);


	template <class T>
	void drop_surface_points_scalar(const UM::Surface& m, UM::PointAttribute<T>& attr, std::string name = "PointAttr", T skip_value = T(-1)) {
		drop_surface(m, name, UM::SurfaceAttributes{ { { "attr", attr.ptr } }, {}, {} }, "attr");
	}

	template <class T>
	void drop_selected_pointset(const UM::PointSet& ps, UM::PointAttribute<T>& attr, std::string name = "PointAttr", T skip_value = T(0)) {
		UM::PointAttribute<bool> sel(ps); for (int v = 0; v < ps.size(); v++)  sel[v] = attr[v] != skip_value;
		UM::Triangles empty_m; empty_m.points = ps;
		drop_surface(empty_m, name, UM::SurfaceAttributes{ { { "selection", sel.ptr } }, {}, {} });
	}

	template <class T>
	void drop_facet_scalar(const UM::Surface& m, UM::FacetAttribute<T>& facet_attr, std::string name = "FacetAttr", T skip_value = T(-1)) {
		if (!drop_mesh_is_active) return;
		UM::Polygons other;
		UM::FacetAttribute<T> other_attr(other);
		other.points = m.points;
		for (int f = 0; f < m.nfacets(); f++) if (facet_attr[f] != skip_value) {
			int nf = other.create_facets(1, m.facet_size(f));
			other_attr[nf] = facet_attr[f];
			for (int fv = 0; fv < m.facet_size(f); fv++)  other.vert(nf, fv) = m.vert(f, fv);
		}
		//drop_surface(m, name, SurfaceAttributes{ {}, { { "attr", facet_attr.ptr } }, {} }, "","attr");
		drop_surface(other, name, UM::SurfaceAttributes{ {}, { { "attr", other_attr.ptr } }, {} }, "", "attr");
	}




	template <class T>
	void drop_corner_scalar(const UM::Surface& m, UM::CornerAttribute<T>& facet_corner_attr, std::string name = "CornerAttr", T skip_value = T(-1)) {
		if (!drop_mesh_is_active) return;
		UM::Triangles outm;
		UM::FacetAttribute<T> attr(outm);
		UM::SurfaceConnectivity fec(m);
		for (int h = 0; h < m.ncorners(); h++) {
			T val = facet_corner_attr[h];
			if (val == skip_value) continue;
			UM::vec3 pts[3] = {
				m.points[m.facets[h]],
				m.points[m.facets[fec.next(h)]],
				m.points[m.facets[fec.prev(h)]]
			};
			int offv = outm.points.create_points(3);
			outm.points[offv] = pts[0];
			outm.points[offv + 1] = 0.7 * pts[0] + .3 * pts[1];
			outm.points[offv + 2] = 0.7 * pts[0] + .3 * pts[2];

			int f = outm.create_facets(1);
			for (int i = 0; i < 3; i++) outm.facets[3 * f + i] = offv + i;
			attr[f] = val;
		}
		drop_surface(outm, name, UM::SurfaceAttributes{ {}, {{ "attr", attr.ptr } }, { } }, "", "attr", "");
	}



	inline UM::Quaternion align_with_uv(UM::vec3 u, UM::vec3 v) {
		v.normalize();
		u.normalize();
		if (u * v < -.99) { UM::Quaternion res;  res.v = UM::vec3(1, 0, 0); res.w = 0;  return res; }
		if (std::abs(u * v) > .99)  return UM::Quaternion();
		UM::vec3 inbetwen_uv(v + u);
		inbetwen_uv.normalize();
		UM::Quaternion res;
		res.w = v * inbetwen_uv; // scalar product with (1,0,0) divided by norm
		res.v = cross(inbetwen_uv, v); // cross product with (1,0,0) 
		return res;
	}

	template <class T>
	void add_arrow(UM::Polygons& out_mesh, UM::FacetAttribute<T>& val, double value, const UM::vec3& A, const UM::vec3& B, ArrowStyle as = ArrowStyle(3, 0)) {
		if (as.diameter == 0) as.diameter = (B - A).norm() / 10.;

		int offv = out_mesh.points.create_points(as.resolution * 3 + 1);
		for (int v = 0; v < as.resolution; v++) out_mesh.points[offv + v] = .5 * as.diameter * UM::vec3(cos(2. * M_PI * double(v) / double(as.resolution)), sin(2. * M_PI * double(v) / double(as.resolution)), 0);
		for (int v = 0; v < as.resolution; v++) out_mesh.points[offv + as.resolution + v] = out_mesh.points[offv + v] + UM::vec3(0, 0, .75);
		for (int v = 0; v < as.resolution; v++) out_mesh.points[offv + 2 * as.resolution + v] = 2. * out_mesh.points[offv + v] + UM::vec3(0, 0, 0.7);
		out_mesh.points[offv + 3 * as.resolution] = UM::vec3(0, 0, 1);
		int offf = out_mesh.create_facets(as.resolution, 4);
		for (int f = 0; f < as.resolution; f++) {
			out_mesh.facets[out_mesh.offset[offf + f] + 0] = offv + f;
			out_mesh.facets[out_mesh.offset[offf + f] + 1] = offv + (f + 1) % as.resolution;
			out_mesh.facets[out_mesh.offset[offf + f] + 2] = offv + (f + 1) % as.resolution + as.resolution;
			out_mesh.facets[out_mesh.offset[offf + f] + 3] = offv + f + as.resolution;
		}
		offf = out_mesh.create_facets(as.resolution, 3);
		for (int f = 0; f < as.resolution; f++) {
			out_mesh.facets[out_mesh.offset[offf + f] + 0] = offv + f + 2 * as.resolution;
			out_mesh.facets[out_mesh.offset[offf + f] + 1] = offv + (f + 1) % as.resolution + 2 * as.resolution;
			out_mesh.facets[out_mesh.offset[offf + f] + 2] = offv + 3 * as.resolution;
		}
		for (int f = 0; f < 2 * as.resolution; f++)  val[offf - as.resolution + f] = T(value);

		UM::vec3 n = B - A;
		double l = std::max(1e-5, n.norm());
		n = n / l;

		for (int v = offv; v < offv + 3 * as.resolution + 1; v++)  out_mesh.points[v][2] *= l;

		UM::Quaternion quat = align_with_uv(UM::vec3(0, 0, 1), n);
		auto M = quat.rotation_matrix();
		//if (n / n.norm() * vec3(0, 0, 1) < -0.99) {
		   // M = mat3x3::identity();
		   // M[2][2] = -1;
		//}
		for (int v = offv; v < offv + 3 * as.resolution + 1; v++)
			out_mesh.points[v] = M * out_mesh.points[v] + A;

	}



	// if as.resolution==0 => draw segments
	// if as.diameter==0 => autoset w.r.t ave_edge_size()
	template <class T>
	void drop_edge_scalar(const UM::Surface& m, UM::CornerAttribute<T>& facet_corner_attr, std::string name = "edge attr", T skip_value = T(-1), ArrowStyle as = ArrowStyle(0, 0)) {
		if (!drop_mesh_is_active) return;
		UM::SurfaceConnectivity fec(m);

		if (as.resolution == 0) {
			UM::PolyLine outm;

			UM::SegmentAttribute<T> attr(outm);
			for (int h = 0; h < m.ncorners(); h++) {
				UM::vec3 G = m.util.bary_verts(fec.facet(h));
				T val = facet_corner_attr[h];
				if (val == skip_value) continue;

				int offv = outm.points.create_points(2);
				outm.points[offv] = 0.9 * m.points[fec.from(h)] + .1 * G;
				outm.points[offv + 1] = 0.9 * m.points[fec.to(h)] + .1 * G;
				int offs = outm.create_segments(1);
				for (int i = 0; i < 2; i++) outm.segments[2 * offs + i] = offv + i;
				attr[offs] = val;
			}
			drop_polyline(outm, name, {}, { { "attr", attr.ptr } }, true);
		}
		else {
			UM::Polygons outm;
			UM::FacetAttribute<T> attr(outm);
			if (as.diameter == 0) as.diameter = .05 * ave_edge_size(m);
			for (int h = 0; h < m.ncorners(); h++) {
				UM::vec3 G = m.util.bary_verts(fec.facet(h));
				T val = facet_corner_attr[h];
				if (val == skip_value) continue;
				add_arrow(outm, attr, facet_corner_attr[h], 0.9 * m.points[fec.from(h)] + .1 * G, 0.9 * m.points[fec.to(h)] + .1 * G, as);
			}
			drop_surface(outm, name, UM::SurfaceAttributes{ {}, { { "attr", attr.ptr } }, {} }, "", "attr");
		}
	}

	// if as.resolution==0 => draw segments
	// if as.diameter==0 => autoset w.r.t ave_edge_size()
	template <class T>
	void drop_edge_scalar_dual(const UM::Surface& m, UM::CornerAttribute<T>& facet_corner_attr, std::string name = "dual edge attr", T skip_value = T(-1), ArrowStyle as = ArrowStyle(0, 0)) {
		if (!drop_mesh_is_active) return;
		UM::SurfaceConnectivity fec(m);
		if (as.resolution == 0) {
			UM::PolyLine outm;
			UM::SegmentAttribute<T> attr(outm);
			for (int h = 0; h < m.ncorners(); h++) {
				UM::vec3 G = m.util.bary_verts(fec.facet(h));

				UM::vec3 Gopp = .5 * (m.points[fec.from(h)] + m.points[fec.to(h)]);
				int opp = fec.opposite(h);
				if (opp != -1) Gopp = m.util.bary_verts(fec.facet(opp));

				T val = facet_corner_attr[h];
				if (val == skip_value) continue;

				int offv = outm.points.create_points(2);
				outm.points[offv] = 0.1 * m.points[fec.from(h)] + .9 * G;
				outm.points[offv + 1] = 0.1 * m.points[fec.from(h)] + .9 * Gopp;
				int offs = outm.create_segments(1);
				for (int i = 0; i < 2; i++) outm.segments[2 * offs + i] = offv + i;
				attr[offs] = val;
			}
			drop_polyline(outm, name, {}, { { "attr", attr.ptr } }, true);
		}
		else {
			UM::Polygons outm;
			UM::FacetAttribute<T> attr(outm);
			if (as.diameter == 0) as.diameter = .05 * ave_edge_size(m);
			for (int h = 0; h < m.ncorners(); h++) {
				UM::vec3 G = m.util.bary_verts(fec.facet(h));

				UM::vec3 Gopp = .5 * (m.points[fec.from(h)] + m.points[fec.to(h)]);
				int opp = fec.opposite(h);
				if (opp != -1) Gopp = m.util.bary_verts(fec.facet(opp));

				T val = facet_corner_attr[h];
				if (val == skip_value) continue;

				add_arrow(outm, attr, facet_corner_attr[h], 0.1 * m.points[fec.from(h)] + .9 * G, 0.1 * m.points[fec.from(h)] + .9 * Gopp, as);
			}

			drop_surface(outm, name, UM::SurfaceAttributes{ {}, { { "attr", attr.ptr } }, {} }, "attr");

		}
	}


	// VOLUMES

	template <class T>
	void extract_surface_scalar(const UM::Volume& m, const UM::CellFacetAttribute<T>& cellfacet_attr, UM::Polygons& surf, UM::FacetAttribute<T>& facet_attr, bool keep_boundary = false, T skip_value = T(-1)) {
		UM::OppositeFacet vc(m);
		for (int c = 0; c < m.ncells(); c++) for (int cf = 0; cf < m.nfacets_per_cell(); cf++) {
			if (!keep_boundary) if (vc.adjacent[m.facet(c, cf)] == -1) continue;
			if (cellfacet_attr[m.facet(c, cf)] == skip_value) continue;
			int off_v = surf.points.create_points(m.facet_size(m.nfacets_per_cell() * c + cf));
			for (int cfv = 0; cfv < m.facet_size(m.nfacets_per_cell() * c + cf); cfv++) surf.points[off_v + cfv] = m.points[m.facet_vert(c, cf, cfv)];
			int off_f = surf.create_facets(1, m.facet_size(m.nfacets_per_cell() * c + cf));
			for (int cfv = 0; cfv < m.facet_size(m.nfacets_per_cell() * c + cf); cfv++) surf.vert(off_f, cfv) = off_v + cfv;
			facet_attr[off_f] = cellfacet_attr[m.facet(c, cf)];

		}
	}




	template <class T>
	void drop_cellfacet_scalar(const UM::Volume& m, UM::CellFacetAttribute<T>& facet_attr, std::string name = "CellFacetAttr", T skip_value = T(-1), bool keep_boundary = false) {
		UM::Polygons surf;
		UM::FacetAttribute<T> surf_attr(surf);
		extract_surface_scalar(m, facet_attr, surf, surf_attr, keep_boundary, skip_value);
		drop_facet_scalar(surf, surf_attr, name, skip_value);
	}

	template <class T, class V = UM::Tetrahedra>
	void drop_volume_points_scalar(const V& m, const UM::PointAttribute<T>& attr, std::string name = "PointAttr", T skip_value = T(-1)) {
		V m_copy;
		m_copy.points.data->assign(m.points.begin(), m.points.end());
		m_copy.cells.assign(m.cells.begin(), m.cells.end());
		UM::PointAttribute<T> attr_copy(m_copy);
		std::vector<bool> tokill(m.nverts(), false);
		for (int v : UM::range(m_copy.nverts())) {
			attr_copy[v] = attr[v];
			tokill[v] = attr[v] == skip_value;
		}
		m_copy.delete_vertices(tokill);
		drop_volume(m_copy, name, UM::VolumeAttributes{ { { "attr", attr_copy.ptr } }, {}, {}, {} }, "", "attr");
	}

	template <class T, class V = UM::Tetrahedra>
	void drop_cells_scalar(const V& m, const UM::CellAttribute<T>& cell_attr, std::string name = "CellAttr", T skip_value = T(-1)) {
		V m_copy;
		m_copy.points.data->assign(m.points.begin(), m.points.end());
		m_copy.cells.assign(m.cells.begin(), m.cells.end());
		UM::CellAttribute<T> attr_copy(m_copy);
		std::vector<bool> tokill(m.ncells(), false);
		for (int c : UM::range(m_copy.ncells())) {
			attr_copy[c] = cell_attr[c];
			tokill[c] = cell_attr[c] == skip_value;
		}
		m_copy.delete_cells(tokill);
		m_copy.delete_isolated_vertices();
		drop_volume(m_copy, name, UM::VolumeAttributes{ {}, { { "attr", attr_copy.ptr } }, {}, {} }, "attr", "");
	}


	template <class T>
	void drop_segments_enlarged(UM::PolyLine& m, const UM::SegmentAttribute<T>& attr, std::string name = "CellAttr", double radius_ratio = 0.1, T skip_value = (T)-1) {
		if (radius_ratio == 0.) {
			drop_polyline(m, name, {}, { { "attr", attr.ptr } }, true);
		}
		else {
			UM::Hexahedra hex;
			hex.create_cells(m.nsegments());
			hex.points.create_points(8 * hex.ncells());
			UM::CellAttribute<T> c_attr(hex);
			for (int c = 0; c < hex.ncells(); c++) {
				c_attr[c] = attr[c];
				UM::vec3 n0 = m.points[m.vert(c, 1)] - m.points[m.vert(c, 0)];
				n0.normalize();
				UM::vec3 n1(n0.y - n0.z, n0.z - n0.x, n0.x - n0.y);
				n1.normalize();
				UM::vec3 n2 = cross(n1, n0);
				n1 = n1 * 0.5 * radius_ratio;
				n2 = n2 * 0.5 * radius_ratio;
				for (int ev = 0; ev < 2; ev++) {
					UM::vec3 v = m.points[m.vert(c, ev)];
					hex.points[8 * c + 4 * ev + 0] = v + n1 + n2;
					hex.points[8 * c + 4 * ev + 1] = v + n1 - n2;
					hex.points[8 * c + 4 * ev + 2] = v - n1 + n2;
					hex.points[8 * c + 4 * ev + 3] = v - n1 - n2;
				}
				for (int j = 0; j < 8; j++) hex.vert(c, j) = 8 * c + j;
			}
			drop_cells_scalar(hex, c_attr, name, skip_value);
		}
	}
};
#endif