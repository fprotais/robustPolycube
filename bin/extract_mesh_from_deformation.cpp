#include <flagging/naive_flagging.h>
#include <utils/trace.h>
#include <filesystem>

#include <deformation/tools.h>
#include <deformation/deformation.h>

#include <utils/LinearSolver.h>

#include "quantization/tools.h"
#include "quantization/blockExtraction.h"
#include "quantization/quantization.h"
#include "quantization/inversion.h"
#include "quantization/chartCleaning.h"

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
using namespace rb_data_structure;
#define FOR(i, n) for(int i = 0; i < n; i++)

int main(int argc, char** argv) {
    Trace::initialize(DEBUG_GRAPHITE_PATH);
    bool has_solver = INITIALISE_LINEARSOLVER(MILP_SOLVER_LPSOLVE, MILP_SOLVER_CPLEX, MILP_SOLVER_GLPK, MILP_SOLVER_GUROBI_CL);
    if (has_solver) std::cerr << "Executable is linked with a MILP solver." << std::endl;
    else std::cerr << "No MILP solver found, will run default algorithm..." << std::endl;

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " mesh.ext flaggingfile polycuboid.ext ElementSizing hexmesh.ext" << std::endl;
        std::cerr << "Input:" << std::endl;
        std::cerr << "- mesh.ext must contain a tet mesh" << std::endl;
        std::cerr << "- flaggingfile must containe a flagging, see README.md and generate_flagging for format info" << std::endl;
        std::cerr << "- polycuboid.ext must contain a polycuboid, same mesh as mesh.ext, with different points" << std::endl;
        std::cerr << "Output: [OPTIONAL]" << std::endl;
        std::cerr << "- ElementSizing, 1 is same hex edge size as polycuboid, 0.5 is twice as small, 2, twice as large, default is 1'" << std::endl;
        std::cerr << "- hexmesh.ext is the extracted hexmesh, default name will be 'mesh_hexmesh.ext'" << std::endl;
        std::cerr << std::endl;
        std::cerr << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cerr << "For more details, see README.md" << std::endl;
        std::cerr << "contact: francoisprotais@gmail.com" << std::endl;
        return 1;
    }
    std::string inputname(argv[1]);
    std::string flaggingname(argv[2]);
    std::string polycuboidname(argv[3]);
    std::string meshname = std::filesystem::path(inputname).stem().string();
    std::string meshextension = std::filesystem::path(inputname).extension().string();
    if (meshextension == ".obj") meshextension = ".mesh";

    std::string hexmeshname = meshname + "_hexmesh" + meshextension;
    double elementsizing = 1;
    if (argc > 4) elementsizing = std::stod(argv[4]);
    if (argc > 5) hexmeshname = argv[5];

    if (elementsizing <= 0) {
        std::cerr << "Elements sizing (=" << elementsizing << ") will lead to degerante elements. We stop here." << std::endl;
        return 1;
    }

    double scaling = 1. / elementsizing;


    // LOADING
    Tetrahedra m;
    Tetrahedra polycuboid;

    // there is some strange manipulation here to account for ultimaille deleting attributes when reading mesh, order of declaration matter.
    read_by_extension(inputname, m);
    read_by_extension(polycuboidname, polycuboid);

    if (m.ncells() == 0) {
        std::cerr << "Mesh given: " << inputname << std::endl << "is empty, we abort." << std::endl;
        return 1;
    }
    if ((m.ncells() != polycuboid.ncells()) || (m.nverts() != polycuboid.nverts())){
        std::cerr << "Meshes given: " << inputname << std::endl << " and "  << polycuboidname << std::endl << "are incompatible, we abort." << std::endl;
        return 1;
    }

    Trace::drop_volume(m, "volume", {});
    Trace::drop_volume(polycuboid, "polycuboid", {});
    
    CellFacetAttribute<int> cfflag(m);

    std::ifstream ifs(flaggingname);
    if (!ifs.is_open()) {
        std::cerr << "Failed opening of flags at : " << flaggingname << std::endl;
        abort();
    }
    FOR(cf, m.ncells() * 4) {
        if (ifs.eof()) cfflag[cf] = -1;
        else ifs >> cfflag[cf];
    }
    ifs.close();

    Trace::drop_cellfacet_scalar(m, cfflag, "flagging", -1, true);

    // RUNNING

    if (!polycuboid_is_valid(m, polycuboid, cfflag)) {
        Trace::alert("The polycuboid has facet not plannar in the flagging??");
        Trace::alert("We are stopping here.");
        Trace::conclude();
        return 1;
    }

    if (!cleanflags(m, polycuboid, cfflag)) {
        Trace::alert("Flagging is locally invalid...");
        Trace::alert("We are stopping here.");
        Trace::drop_cellfacet_scalar(m, cfflag, "finalFlaggingWithTries", -1, true);
        Trace::conclude();
        return 1;
    }
    Trace::drop_cellfacet_scalar(m, cfflag, "corrected_flagging", -1, true);


    Sorted_charts charts;

    compute_charts(m, polycuboid, cfflag, charts);
    disp_charts_dim(m, charts, 0, "charts");
    disp_charts_dim(m, charts, 1, "charts");
    disp_charts_dim(m, charts, 2, "charts");
    Block_decomposition blocks;
    chart_voxelisation(m, polycuboid, charts, blocks);

    quantize(charts, blocks, scaling);

    disp_block_decomposition(blocks, "Blocks");
    Final_mesh hexmesh;

    inverse(m, polycuboid, charts, blocks, hexmesh);

    Trace::drop_cells_scalar(hexmesh.m, hexmesh.original_bloc, "hexmesh");
    Trace::drop_cellfacet_scalar(hexmesh.m, hexmesh.bnd_chart, "hexmesh_charts", -1, true);
    // SAVING
    write_by_extension(hexmeshname, hexmesh.m, { {},{},{},{} });

    Trace::conclude();
    return 0;
}