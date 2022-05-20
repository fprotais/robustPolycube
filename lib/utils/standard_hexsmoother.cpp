#include "standard_hexsmoother.h"


using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

bool smooth_tet_mesh(UM::Tetrahedra& m, std::vector<bool>& locks, const smoother_options& options) {
    assert(m.nverts() * 3 == locks.size());

    std::vector<double> verts(m.nverts() * 3);
    FOR(v, m.nverts()) FOR(d, 3) verts[3 * v + d] = m.points[v][d];
    std::vector<std::array<int, 4>> tets(m.ncells());
    FOR(f, m.ncells()) FOR(fv, 4) tets[f][fv] = m.vert(f, fv);
    Tets_id_with_lock var(verts, tets, locks);

    Elliptic_smoother_3D opt(var, tets.size(), options);
    bool res = opt.go();

    var.get_verts(verts);
    FOR(v, m.nverts()) FOR(d, 3) m.points[v][d] = verts[3 * v + d];
    return res;
}
bool smooth_tet_mesh(UM::Tetrahedra& m, const UM::Tetrahedra& ref, std::vector<bool>& locks, const smoother_options& options){
    assert(m.nverts() == ref.nverts());
    assert(m.nverts() * 3 == locks.size());

    std::vector<double> verts(m.nverts() * 3);
    FOR(v, m.nverts()) FOR(d, 3) verts[3 * v + d] = m.points[v][d];
    std::vector<std::array<int, 4>> tets(m.ncells());
    FOR(f, m.ncells()) FOR(fv, 4) tets[f][fv] = m.vert(f, fv);
    Tets_id_with_lock var(verts, tets, locks);

    Elliptic_smoother_3D opt(var, tets.size(), options);
    FOR(c, tets.size()) {
        std::array<UM::vec3, 4> tet_ref;
        FOR(cv, 4) tet_ref[cv] = ref.points[tets[c][cv]];
        opt.set_tet_ref(c, tet_ref);
    }
    bool res = opt.go();

    var.get_verts(verts);
    FOR(v, m.nverts()) FOR(d, 3) m.points[v][d] = verts[3 * v + d];
    return res;
}

constexpr int HEX_CORNER_SPLITTING[8][4] = {
	{0,1,2,4}, {1,3,0,5}, {2,0,3,6}, {3,2,1,7},
	{4,6,5,0}, {5,4,7,1}, {6,7,4,2}, {7,5,6,3},
};

constexpr int TET_MATRIX_HEX_LEXICOGRAPHIC_SPLIT[24][4] = {
    {0,1,3,4},{0,1,2,5},{0,2,6,1},{0,2,4,3},
    {0,4,1,6},{0,4,5,2},{1,3,0,7},{1,3,2,5},
    {1,5,3,4},{1,5,7,0},{2,3,7,0},{2,3,6,1},
    {2,6,0,7},{2,6,4,3},{3,7,6,1},{3,7,2,5},
    {4,5,0,7},{4,5,1,6},{4,6,5,2},{4,6,7,0},
    {5,7,1,6},{5,7,3,4},{6,7,5,2},{6,7,4,3}
};



bool smooth_hex_mesh(UM::Hexahedra& m, const std::vector<bool>& locks) {
	std::cerr << "Running hex smoother...";
	assert(m.nverts() * 3 == locks.size());
	smoother_options options = _3D_default;
	options.static_threshold = 1e-4;
	options.bfgs_threshold = 1e-8;
	options.theta = 1e-3;
	options.bfgs_maxiter = 100;
	options.eps_from_theorem = false;
	options.maxiter = 10;
	options.debug = false;
	std::vector<double> verts(m.nverts() * 3);
	FOR(v, m.nverts()) FOR(d, 3) verts[3 * v + d] = m.points[v][d];

	Tetrahedra proxy_mesh;
	proxy_mesh.points = m.points;
	std::vector<std::array<UM::vec3, 4>> refs;

	std::vector<int> map(m.nverts(), -1);
	const std::array<vec3, 8> default_cell = { { {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {1.,1.,0.}, {0.,0.,1.}, {1.,0.,1.}, {0.,1.,1.}, {1.,1.,1.} } };
	FOR(c, m.ncells()) {
		std::array<bool, 8> has_on_in = { 0 };
		int start_c = proxy_mesh.ncells();
		FOR(i, 8) FOR(j, 4) FOR(d, 3) if (!locks[3 * m.vert(c, HEX_CORNER_SPLITTING[i][j]) + d]) has_on_in[i] = true;
		FOR(i, 8) if (has_on_in[i]) {
			int new_c = proxy_mesh.create_cells(1);
			FOR(j, 4) proxy_mesh.vert(new_c, j) = m.vert(c, HEX_CORNER_SPLITTING[i][j]);
		}
		FOR(i, 8) map[m.vert(c, i)] = i;
		int nb_ = 0;
		FOR(i, 8) if (has_on_in[i]) {
			std::array<UM::vec3, 4> tet_ref;
			FOR(cv, 4) tet_ref[cv] = default_cell[map[proxy_mesh.vert(start_c + nb_, cv)]];
			refs.push_back(tet_ref);
			nb_++;
		}
	}



	std::vector<std::array<int, 4>> tets(proxy_mesh.ncells());
	FOR(f, proxy_mesh.ncells()) FOR(fv, 4) tets[f][fv] = proxy_mesh.vert(f, fv);
	Tets_id_with_lock var(verts, tets, locks);

	Elliptic_smoother_3D opt(var, tets.size(), options);
	FOR(i, refs.size()) opt.set_tet_ref(i, refs[i]);

	bool res = opt.go();
	var.get_verts(verts);
	FOR(v, m.nverts()) FOR(d, 3) m.points[v][d] = verts[3 * v + d];
	std::cerr << "Done\n";
	return res;
}
