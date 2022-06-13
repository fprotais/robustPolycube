#include <utils/trace.h>
#include <filesystem>


#include "post_processing/pillowing.h"
#include "post_processing/CADaware_hexsmoothing.h"


#ifndef DEBUG_GRAPHITE_PATH
#define DEBUG_GRAPHITE_PATH "C:/fprotais/softwares/graphite/build/Windows/bin/Release/graphite.exe"
#endif 

#ifndef MILP_SOLVER_CPLEX
#define MILP_SOLVER_CPLEX "cplex"
#endif 
#ifndef MILP_SOLVER_LPSOLVE
#define MILP_SOLVER_LPSOLVE ""
#endif 
#ifndef MILP_SOLVER_GUROBI_CL
#define MILP_SOLVER_GUROBI_CL ""
#endif 
#ifndef MILP_SOLVER_GLPK
#define MILP_SOLVER_GLPK ""
#endif 

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

int main(int argc, char** argv) {
    Trace::initialize(DEBUG_GRAPHITE_PATH);

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " mesh.ext hexmesh.ext result.ext" << std::endl;
        std::cerr << "Input:" << std::endl;
        std::cerr << "- mesh.ext must contain a tet mesh" << std::endl;
        std::cerr << "- hexmesh.ext must containe a hex mesh" << std::endl;
        std::cerr << "-> Boundary of hexmesh and mesh must be pretty similar" << std::endl;
        std::cerr << "Output: [OPTIONAL]" << std::endl;
        std::cerr << "- result.ext is the smoothed hexmesh, default name will be 'hexmesh_smoothed.ext'" << std::endl;
        std::cerr << std::endl;
        std::cerr << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cerr << "For more details, see README.md" << std::endl;
        std::cerr << "contact: francoisprotais@gmail.com" << std::endl;
        return 1;
    }
    std::string inputname(argv[1]);
    std::string hexmeshname(argv[2]);
    std::string meshname = std::filesystem::path(inputname).stem().string();
    std::string meshextension = std::filesystem::path(inputname).extension().string();

    std::string resultname = meshname + "_smoothed" + meshextension;
    if (argc > 3) resultname = argv[3];


    // LOADING
    Tetrahedra m;
    Hexahedra hexmesh;

    // there is some strange manipulation here to account for ultimaille deleting attributes when reading mesh, order of declaration matter.
    read_by_extension(inputname, m);
    read_by_extension(hexmeshname, hexmesh);

    if (m.ncells() == 0) {
        std::cerr << "Mesh given: " << inputname << std::endl << "is empty, we abort." << std::endl;
        return 1;
    }
    if (hexmesh.ncells() == 0) {
        std::cerr << "Mesh given: " << hexmeshname << std::endl << "is empty, we abort." << std::endl;
        return 1;
    }


    Trace::drop_volume(m, "volume", {});
    Trace::drop_volume(hexmesh, "input_hexmesh", {});
    
    Triangles bnd;
    OppositeFacet vec(m);
    bnd.points.data->assign(m.points.begin(), m.points.end());
    FOR(cf, 4 * m.ncells()) if (vec.adjacent[cf] == -1) {
        int f = bnd.create_facets(1);
        FOR(fv, 3) bnd.vert(f, fv) = m.facet_vert(cf / 4, cf % 4, fv);
    }
    bnd.delete_isolated_vertices();
    Trace::drop_surface(bnd, "surface", {});

    // RUNNING
    CellAttribute<int> blocks_id(hexmesh);
    add_a_pillow(hexmesh, blocks_id);
    Trace::drop_cells_scalar(hexmesh, blocks_id, "pillowed", 1);

    smooth(hexmesh, bnd, 1);
    Trace::drop_cells_scalar(hexmesh, blocks_id, "smoothed", 1);

    // SAVING
    write_by_extension(resultname, hexmesh, { {},{},{},{} });

    Trace::conclude();
    return 0;
}