#include "trace.h"
#include <iostream>
#include <fstream>
#include <map>

#define FOR(i, n) for(int i = 0; i < n; i++)
using namespace UM;

struct Args {
	static bool has(const std::string& s) {
		return data.find(s) != data.end();
	}

	static void add_param(const std::string& s) {
		int eq_pos = 0;
		FOR(i, s.size()) if (s[i] == '=') eq_pos = i;
		data[s.substr(0, eq_pos)] = s.substr(eq_pos + 1, s.size());
	}

	static void init(int argc, char** argv) {
		FOR(p, argc - 1) add_param(argv[p + 1]);
		if (data.size() == 1 && data.find("args_file") != data.end())  init(data["args_file"]);
		std::cerr << "Binary " << argv[0] << " called with Args:\n ";
		show();
	}

	static void init(const std::string& filename) {
		std::ifstream infile(filename);
		std::string s;
		while (std::getline(infile, s))  add_param(s);
	}

	static void show() {
		std::cerr << "------------------Global parameters------------------\n";
		for (auto it = data.begin(); it != data.end(); it++)
			std::cerr << it->first << " => " << it->second << std::endl;
	}

	template <class T> static T get(const std::string& name) {
		if constexpr (std::is_same_v<T, std::string>)  return data[name];
		if constexpr (std::is_same_v<T, bool>)         return data[name] == "True";
		if constexpr (std::is_same_v<T, int>)          return std::stoi(data[name]);
		if constexpr (std::is_same_v<T, unsigned int>) return std::stoul(data[name]);
		if constexpr (std::is_same_v<T, float>)        return std::stof(data[name]);
		if constexpr (std::is_same_v<T, double>)       return std::stod(data[name]);
		um_assert(false);
	}

	static std::map<std::string, std::string> data;
};

inline double det(vec2 a, vec2 b) { return a.x * b.y - b.x * a.y; }

inline double vector_angle(vec3 v0, vec3 v1) { return atan2(cross(v0, v1).norm(), v0 * v1); }
inline double corner_angle(const SurfaceConnectivity& fec, int c) { return vector_angle(fec.geom(c), -fec.geom(fec.prev(c))); }

double ave_edge_size(const Surface& m) {
	SurfaceConnectivity fec(m);
	double res = 0;
	FOR(c, m.ncorners())  res += fec.geom(c).norm();
	return res / double(m.ncorners());
}



#ifdef WIN32
#include <windows.h>
void color(int t, int f) {
	HANDLE H = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(H, f * 16 + t);
}
#else
void color(int t, int f) {
	if (t == 15 && f == 0) {
		std::cerr << "\e[0m";
		return;
	}
	std::cerr << "\e[0";
	if (t & 8) std::cerr << ";1";
	t = ((t & 1) << 2) | (t & 2) | ((t >> 2) & 1);
	std::cerr << ";3" << char('0' + char(t));
	f = ((f & 1) << 2) | (f & 2) | ((f >> 2) & 1);
	if (f != 0) std::cerr << ";4" << char('0' + char(f));
	std::cerr << 'm';
}
#endif
using namespace UM;

namespace Trace {

	bool drop_mesh_is_active = false;
	static bool trace_steps_active = true;
	static int num_drop = 0;
	static std::string outputdir = "";
	static std::string outputprefix = "";
	static std::string graphite_path;


	struct  LogTime {
		/**
		 * CheckPoint is the list item of LogTime
		 * "up" is it's father
		 * "right" is the next item at the same level
		 */
		struct CheckPoint {
			CheckPoint(std::string const& p_n, unsigned int p_up) { n = p_n; up = p_up; t = clock(); right = (unsigned int)(-1); }
			std::string n;
			clock_t t;
			unsigned int right;
			unsigned int up;
		};

		static std::vector<CheckPoint> check;
		static std::vector<std::pair<std::string, double> > out_values;
		static std::vector<std::pair<std::string, std::string> > out_strings;


		static bool is_start_section(unsigned int i) {
			return check[i].right != i + 1;
		}

		static bool is_end_section(unsigned int i) {
			return check[i].n == "end section";
		}
		static bool is_final(unsigned int i) {
			return i + 1 == check.size();
		}

		static double time(unsigned int i) {
			return (double(check[check[i].right].t) - double(check[i].t)) / double(CLOCKS_PER_SEC);
		}

		static unsigned int dec(unsigned int i = (unsigned int)(-1)) {
			if (i == (unsigned int)(-1))
				i = (unsigned int)(check.size() - 1);
			unsigned int res = 0;
			i = check[i].up;
			while (i != (unsigned int)(-1)) {
				res++;
				i = check[i].up;
			}
			return res;
		}

		static unsigned int lastdec() {
			if (!check.empty())
				return dec((unsigned int)(check.size() - 1));
			return 0;
		}

		static void debug() {
			std::cerr << std::endl;
			std::cerr << "--------------BEGIN DEBUG-------------------" << std::endl;
			for (size_t i = 0; i < check.size(); i++) {
				std::cerr << std::string((unsigned int)(4 * dec((unsigned int)(i))), ' ') << i << "  r = " << check[i].right << " u = " << check[i].up
					<< "\tstart" << check[i].t << "\tname" << check[i].n << std::endl;
			}
			std::cerr << std::endl;
			std::cerr << "-------------- END  DEBUG-------------------" << std::endl;
		}

		static std::string cur_stack() {
			std::string res;
			std::vector<unsigned int> stack;
			{
				unsigned int i = (unsigned int)(check.size() - 1);
				while (i != (unsigned int)(-1)) { stack.push_back(i); i = check[i].up; }
			}
			for (int i = int(stack.size()) - 1; i >= 0; i--) {
				res.append(check[stack[size_t(i)]].n);
				if (i > 0) res.append(" ==> ");
			}
			return res;
		}

		static void report(std::ostream& out, unsigned int timing_depth = 10000) {
			if (check.empty()) return;
			if (check.back().n != "the end") { add_step("the end"); }

			if (timing_depth != (unsigned int)(-1)) {
				out << "\n***********************************************************" << std::endl;
				out << "                  TIMING SUMMARY " << std::endl;
				for (unsigned int i = 0; i < check.size() - 1; i++) {
					if (dec(i) > timing_depth) continue;
					if (is_start_section(i))
						out << std::string(4U * dec(i), ' ') << time(i) << "\t====>  " << check[i].n << std::endl;
					else if (!is_end_section(i) && check[i].n != "begin section")
						out << std::string(4U * dec(i), ' ') << time(i) << "\t" << check[i].n << std::endl;
				}

				out << double(check[check.size() - 1U].t - check[0].t) / double(CLOCKS_PER_SEC) << "\tTOTAL" << std::endl;
			}
			out << "\n***********************************************************" << std::endl;
			out << "                  OUPUT VALUES" << std::endl;
			for (auto const& v : out_values)
				out << v.second << " \t" << v.first << std::endl;
		}



		static void report_py(std::ostream& out, unsigned int timing_depth) {
			if (!check.empty() && check.back().n != "the end") { add_step("the end"); }


			out << "{\"finished\": \"yes\"";
			if (timing_depth != (unsigned int)(-1)) {
				for (unsigned int i = 0; i < check.size() - 1; i++) {
					if (dec(i) > timing_depth) continue;
					if (is_start_section(i)) out << ",\"TIME_" << check[i].n << "\" :  " << time(i);
					else if (!is_end_section(i) && check[i].n != "begin section")
						out << ",\"TIME_" << check[i].n << "\":  " << time(i);
				}
			}
			for (auto const& v : out_values)
				out << ", \"" << v.first << "\": " << v.second;
			for (auto const& v : out_strings)
				out << ",\"" << v.first << "\": \"" << v.second << "\"";
			out << "}";

		}





		// construct API
		static void log_value(std::string const& str, double val) {
			std::cerr << "LogTime >> " << str << " = " << val << std::endl;
			out_values.push_back(std::pair<std::string, double>(str, val));
		}
		static void log_string(std::string const& str, std::string const& val) {
			std::cerr << "LogTime >> " << str << " = " << val << std::endl;
			out_strings.push_back(std::pair<std::string, std::string>(str, val));
		}

		static void add_step(std::string const& name) {
			if (check.empty()) {
				CheckPoint c(name, (unsigned int)(-1));
				check.push_back(c);
			}
			else {
				CheckPoint c(name, check.back().up);
				check.back().right = (unsigned int)(check.size());
				check.push_back(c);
			}
			const char* symbol = "#>=_-~..............";
			std::cerr << std::endl << std::string(4 * lastdec(), ' ') << std::string(80 - 4 * lastdec(), symbol[dec()]) << std::endl;
			std::cerr << std::string(4 * lastdec() + 4, ' ') << cur_stack() << std::endl << std::endl;
		}

		static void start_section(std::string const& secname, std::string const& name = "begin section") {
			add_step(secname);
			CheckPoint c(name, (unsigned int)(check.size()) - 1);
			check.push_back(c);
		}
		static void end_section() {
			unsigned int u = check.back().up;
			check[u].right = (unsigned int)(check.size());
			check.back().right = (unsigned int)(check.size());
			CheckPoint c("end section", check[u].up);
			check.push_back(c);
		}


		static void drop_file(std::string const& filename, bool append, unsigned int timing_depth) {
			std::ofstream f;
			if (append) f.open(filename.c_str(), std::fstream::app);
			else		f.open(filename.c_str());
			report(f, timing_depth);
			f.close();
		}
	};

	std::vector<LogTime::CheckPoint> LogTime::check;
	std::vector<std::pair<std::string, double> > LogTime::out_values;
	std::vector<std::pair<std::string, std::string> > LogTime::out_strings;








	void SwitchDropInScope::force_activity(bool b) { drop_mesh_is_active = b; }

	SwitchDropInScope::SwitchDropInScope(bool make_active) {
		if (outputdir == "") make_active = false;
		save_val = drop_mesh_is_active;
		drop_mesh_is_active = make_active;
	};
	SwitchDropInScope::~SwitchDropInScope() {
		drop_mesh_is_active = save_val;
	}

	void SwitchTextInScope::force_activity(bool b) { trace_steps_active = b; }
	SwitchTextInScope::SwitchTextInScope(bool make_active) {
		if (outputdir == "") make_active = false;
		save_val = trace_steps_active;
		trace_steps_active = make_active;
	};
	SwitchTextInScope::~SwitchTextInScope() {
		trace_steps_active = save_val;
	}

	void initialize(const std::string& graphite_path_in) {
		outputdir = "";
		outputprefix = "debug_";
		if (graphite_path == "") {
			std::remove("view.lua");
			drop_mesh_is_active = true;
			trace_steps_active = true;
			graphite_path = graphite_path_in;
		}
	};
	void conclude() {
		if (drop_mesh_is_active) {
			std::string cmd = graphite_path + " " + outputdir + "view.lua";
			system(cmd.c_str());
		}
	}


	void set_alert_level(int alert_level) {
		if (alert_level == 0) color(10, 0);
		if (alert_level == 1) color(14, 0);
		if (alert_level == 2) color(12, 0);
	}
	void step(std::string stepname, int alert_level) {
		if (!trace_steps_active) return;
		set_alert_level(alert_level);
		LogTime::add_step(stepname);
		color(15, 0);
	}
	Section::Section(const std::string& str) {
		trace_was_active = Trace::trace_steps_active;
		if (!trace_was_active) return;
		set_alert_level(0);
		LogTime::start_section(str);
		color(15, 0);
	}
	Section::~Section() { if (trace_was_active) LogTime::end_section(); }

	void log_value(std::string const& str, double val, int alert_level) { set_alert_level(alert_level); LogTime::log_value(str, val); color(15, 0); }
	void log_string(std::string const& str, std::string const& val, int alert_level) { set_alert_level(alert_level); LogTime::log_string(str, val); color(15, 0); }
	void alert(std::string msg) { log_string("ALERT", msg, 2); }
	void show_log() {
		LogTime::report(std::cerr);
	}
	void append_py_log(std::string const& filename) {
		std::ofstream myfile;
		//myfile.open(filename, std::ios::out | std::ios::app);
		myfile.open(filename, std::ios::out);
		LogTime::report_py(myfile, -1);
		myfile.close();
	}


	void drop_volume(const Volume& m,
		const std::string& name,
		const VolumeAttributes& attr,
		std::string cell_attr,
		std::string point_attr
	) {
		if (!drop_mesh_is_active) return;
		std::string filename = outputprefix + name + "_" + std::to_string(num_drop++) + ".geogram";
		write_by_extension(filename, m, attr);
		std::ofstream myfile;
		myfile.open(outputdir + "view.lua", std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		if (!cell_attr.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'cells." << cell_attr << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
		}
		else if (!point_attr.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'vertices." << point_attr << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
			myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 3'\n";
		}
		myfile.close();
	}


	void drop_surface(const Surface& m,
		const std::string& name,
		const SurfaceAttributes& attr,
		std::string point_attr,
		std::string facet_attr,
		std::string corner_attr
	) {
		if (!drop_mesh_is_active) return;
		std::string filename = outputprefix + name + "_" + std::to_string(num_drop++) + ".geogram";
		write_by_extension(filename, m, attr);
		std::ofstream myfile;
		myfile.open(outputdir + "view.lua", std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		myfile << "scene_graph.current().shader.mesh_style = 'true; 0 0 0 1; 1'\n";

		if (!corner_attr.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'facet_corners." << corner_attr << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
		}
		else if (!facet_attr.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'facets." << facet_attr << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
		}
		else if (!point_attr.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'vertices." << point_attr << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
			myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 3'\n";
		}
		else {
			myfile << "scene_graph.current().shader.mesh_style = 'true; 0 0 0 1; 1'\n";

		}
		myfile.close();
	}


	void drop_texture(const Surface& m, const std::string& name, CornerAttribute<vec2>& u) {
		if (!drop_mesh_is_active) return;
		std::string filename = outputprefix + name + "_" + std::to_string(num_drop++) + ".geogram";
		write_by_extension(filename, m, SurfaceAttributes{ {}, {}, { { "U", u.ptr } } });
		std::ofstream myfile;
		myfile.open(outputdir + "view.lua", std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		myfile << "scene_graph.current().shader.painting = 'TEXTURE'\n";
		myfile << "scene_graph.current().shader.tex_coords = 'facet_corners.U'\n";
		myfile.close();
	}



	void drop_polyline(const PolyLine& m,
		const std::string& name,
		std::vector<NamedContainer> vertex_attributes,
		std::vector<NamedContainer> segment_attributes,
		bool show_edge_attr) {
		if (!drop_mesh_is_active) return;
		std::string filename = outputprefix + name + "_" + std::to_string(num_drop++) + ".geogram";
		write_by_extension(filename, m, PolyLineAttributes{ vertex_attributes, segment_attributes });
		std::ofstream myfile;
		myfile.open(outputdir + "view.lua", std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		if (show_edge_attr) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'edges.attr'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
		}
		myfile.close();
	}


	void drop_facet_vec3(const Surface& m, FacetAttribute<vec3>& facet_attr, std::string name, ArrowStyle as, bool proportional_to_ave_edge_size, bool normalize) {
		double scale = 1;
		if (proportional_to_ave_edge_size) scale = ave_edge_size(m);
		if (as.resolution == 0) {
			PolyLine outm;

			FOR(f, m.nfacets()) {
				vec3 G = m.util.bary_verts(f);
				int offv = outm.points.create_points(2);
				outm.points[offv] = G;
				vec3 decal = facet_attr[f];
				if (normalize) decal.normalize();
				decal = scale * decal;
				outm.points[offv + 1] = G + decal;
				int offs = outm.create_segments(1);
				FOR(i, 2) outm.segments[2 * offs + i] = offv + i;
			}
			drop_polyline(outm, name, {}, {});
		}
		else {
			Polygons outm;
			FacetAttribute<double> attr(outm);
			if (as.diameter == 0) as.diameter = .05 * ave_edge_size(m);
			FOR(f, m.nfacets()) {
				vec3 G = m.util.bary_verts(f);
				vec3 decal = facet_attr[f];
				if (normalize) decal.normalize();
				decal = scale * decal;
				add_arrow(outm, attr, facet_attr[f].norm(), G, G + decal, as);
			}
			drop_surface(outm, name, SurfaceAttributes{ {}, { { "attr", attr.ptr } }, {} }, "", "attr");
		}
	}
	void drop_surface_points_vec3(const Surface& m, PointAttribute<vec3>& attr, std::string name, ArrowStyle as, bool proportional_to_ave_edge_size, bool normalize) {
		double scale = 1;
		if (proportional_to_ave_edge_size) scale = ave_edge_size(m);
		if (as.resolution == 0) {
			PolyLine outm;
			FOR(v, m.nverts()) {
				vec3 G = m.points[v];
				int offv = outm.points.create_points(2);
				outm.points[offv] = G;
				vec3 decal = attr[v];
				if (normalize) decal.normalize();
				decal = scale * decal;
				outm.points[offv + 1] = G + decal;
				int offs = outm.create_segments(1);
				FOR(i, 2) outm.segments[2 * offs + i] = offv + i;
			}
			drop_polyline(outm, name, {}, {});
		}
		else {
			Polygons outm;
			FacetAttribute<double> attrf(outm);
			if (as.diameter == 0) as.diameter = .05 * ave_edge_size(m);
			FOR(v, m.nverts()) {
				vec3 G = m.points[v];
				vec3 decal = attr[v];
				if (decal.norm2() < 1e-20) continue;
				if (normalize) decal.normalize();
				decal = scale * decal;
				add_arrow(outm, attrf, attr[v].norm(), G, G + decal, as);
			}

			drop_surface(outm, name, SurfaceAttributes{ {}, { { "attrf", attrf.ptr } }, {} }, "", "attr");
		}
	}


};