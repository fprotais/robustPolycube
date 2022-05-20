#include <utils/trace.h>
#include <filesystem>

#include <flagging/naive_flagging.h>

#include <deformation/tools.h>
#include <deformation/deformation.h>

#include <utils/LinearSolver.h>

#include "quantization/tools.h"
#include "quantization/blockExtraction.h"
#include "quantization/quantization.h"
#include "quantization/inversion.h"

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
using namespace rb_data_structure;
#define FOR(i, n) for(int i = 0; i < n; i++)

int main(int argc, char** argv) {
    Trace::initialize(DEBUG_GRAPHITE_PATH);
    bool has_solver = INITIALISE_LINEARSOLVER(MILP_SOLVER_LPSOLVE, MILP_SOLVER_CPLEX, MILP_SOLVER_GLPK, MILP_SOLVER_GUROBI_CL);
    if (has_solver) std::cerr << "Executable is linked with a MILP solver." << std::endl;
    else std::cerr << "No MILP solver found, will run default algorithm..." << std::endl;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " mesh.ext ElementSizing hexmesh.ext flaggingfile" << std::endl;
        std::cerr << "Input:" << std::endl;
        std::cerr << "- mesh.ext must contain a tet or tri mesh" << std::endl;
        std::cerr << "Output: [OPTIONAL]" << std::endl;
        std::cerr << "- ElementSizing, 1 is same hex edge size as mesh, 0.5 is twice as small, 2, twice as large, default is 1'" << std::endl;
        std::cerr << "- hexmesh.ext is the extracted hexmesh, default name will be 'mesh_hexmesh.ext'" << std::endl;
        std::cerr << "- flaggingfile is a file containing a flagging, without it, we compute a naive flagging" << std::endl;
        std::cerr << std::endl;
        std::cerr << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cerr << "For more details, see README.md" << std::endl;
        std::cerr << "contact: francoisprotais@gmail.com" << std::endl;
        return 1;
    }
    std::string inputname(argv[1]);
    std::string meshname = std::filesystem::path(inputname).stem().string();
    std::string meshextension = std::filesystem::path(inputname).extension().string();
    if (meshextension == ".obj") meshextension = ".mesh";

    std::string hexmeshname = meshname + "_hexmesh" + meshextension;
    std::string flaggingname = "";

    double elementsizing = 1;
    if (argc > 2) elementsizing = std::stod(argv[2]);
    if (argc > 3) hexmeshname = argv[3];
    if (argc > 4) flaggingname = argv[4];
    

    if (elementsizing <= 0) {
        std::cerr << "Elements sizing (=" << elementsizing << ") will lead to degerante elements. We stop here." << std::endl;
        return 1;
    }

    double scaling = 1. / elementsizing;


    // LOADING and HANDLING FLAGGING
    Tetrahedra volume;
    Triangles surf;

    // there is some strange manipulation here to account for ultimaille deleting attributes when reading mesh, order of declaration matter.
    read_by_extension(inputname, volume);
    read_by_extension(inputname, surf);
    FacetAttribute<int> flag(surf);

    if (volume.ncells() != 0) {
        Trace::drop_volume(volume, "volume", {});
        CellFacetAttribute<int> cflag(volume, -1);

        if (flaggingname != "") {
            std::ifstream ifs(flaggingname);
            if (!ifs.is_open()) {
                std::cerr << "Failed opening of flags at : " << flaggingname << std::endl;
                abort();
            }
            FOR(cf, volume.ncells() * 4) {
                if (ifs.eof()) cflag[cf] = -1;
                else ifs >> cflag[cf];
            }
            ifs.close();
        }
        else {
            generate_naive_flagging(volume, cflag);
        }

        transfer_to_surf(volume, cflag, surf, flag);
    }
    else {
        if (surf.nfacets() == 0) {
            std::cerr << "Mesh given: " << inputname << std::endl << "is empty, we cannot compute deformation." << std::endl;
            return 1;
        }
        Trace::drop_surface(surf, "surface", {});
        if (flaggingname != "") {

            std::ifstream ifs(flaggingname);
            if (!ifs.is_open()) {
                std::cerr << "Failed opening of flags at : " << flaggingname << std::endl;
                abort();
            }
            FOR(f, surf.nfacets()) {
                if (ifs.eof()) flag[f] = -1;
                else ifs >> flag[f];
            }
            ifs.close();
        }
        else {
            generate_naive_flagging(surf, flag);
        }
    }


    //  DEFORMATION

    center_and_normalise_mesh(surf);
    Trace::drop_facet_scalar(surf, flag, "flagging");

    Tetrahedra m;
    CellAttribute<int> inside(m);
    CellFacetAttribute<int> cfflags(m);

    put_model_in_box(surf, flag, m, cfflags, inside);

    Trace::drop_cellfacet_scalar(m, cfflags, "volume_flagging", -1, true);
    Trace::drop_cells_scalar(m, inside, "embedded_mesh");
    PointAttribute<vec3> U(m);
    cube_cover(m, cfflags, U);
    correct_param(m, cfflags, U);

    remove_outerbox(m, cfflags, inside);
    disp_polycube(m, U, cfflags, "corrected");

    Tetrahedra polycuboid;
    *polycuboid.points.data = U.ptr->data;
    polycuboid.cells = m.cells;
    
    // QUANTIZATION

    cleanflags(m, polycuboid, cfflags);

    Sorted_charts charts;

    compute_charts(m, polycuboid, cfflags, charts);
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

    // POST PRoCESSING

    add_a_pillow(hexmesh.m, hexmesh.original_bloc);
    Trace::drop_cells_scalar(hexmesh.m, hexmesh.original_bloc, "pillowed", -2);

    smooth(hexmesh.m, surf, 1);
    Trace::drop_cells_scalar(hexmesh.m, hexmesh.original_bloc, "smoothed", 1);


    // SAVING
    write_by_extension(hexmeshname, hexmesh.m, { {},{},{},{} });

    Trace::conclude();
    return 0;
}