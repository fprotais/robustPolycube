#include <utils/trace.h>
#include <filesystem>

#include <deformation/tools.h>
#include <deformation/deformation.h>

#ifndef DEBUG_GRAPHITE_PATH
#define DEBUG_GRAPHITE_PATH "C:/fprotais/softwares/graphite/build/Windows/bin/Release/graphite.exe"
#endif 

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

int main(int argc, char** argv) {
    Trace::initialize(DEBUG_GRAPHITE_PATH);
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " mesh.ext flaggingfile remeshed.ext remeshedflagging polycuboid.ext" << std::endl;
        std::cerr << "Input:" << std::endl;
        std::cerr << "- mesh.ext must contain a tet mesh or tri mesh" << std::endl;
        std::cerr << "- flaggingfile must contain a flagging, see README.md and generate_flagging for format info" << std::endl;
        std::cerr << "Output: [OPTIONAL]" << std::endl;
        std::cerr << "- remeshed.ext is the inner mesh of mesh.ext, if a tet mesh is given as input, /!\\ we remesh it /!\\, see README.md for details, default name will be 'mesh_remeshed.ext'" << std::endl;
        std::cerr << "- remeshedflagging flagging on remeshed, default name will be 'mesh_remeshed.flags'" << std::endl;
        std::cerr << "- polycuboid.ext is result deformed mesh -polycuboid- of remeshed.ext, see README.md for details, default name will be 'mesh_polycuboid.ext'" << std::endl;
        std::cerr << std::endl;
        std::cerr << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cerr << "For more details, see README.md" << std::endl;
        std::cerr << "contact: francoisprotais@gmail.com" << std::endl;
        return 1;
    }
    std::string inputname(argv[1]);
    std::string flaggingname(argv[2]);
    std::string meshname = std::filesystem::path(inputname).stem().string();
    std::string meshextension = std::filesystem::path(inputname).extension().string();
    if (meshextension == ".obj") meshextension = ".mesh";

    std::string remeshedname = meshname + "_remeshed" + meshextension;
    std::string remeshedflaggingname = meshname + "_remeshed.flags";
    std::string polycuboidname = meshname + "_polycuboid" + meshextension;
    if (argc > 3) remeshedname = argv[3];
    if (argc > 4) remeshedflaggingname = argv[4];
    if (argc > 5) polycuboidname = argv[5];

    // LOADING

    Tetrahedra volume;
    Triangles surf;

    // there is some strange manipulation here to account for ultimaille deleting attributes when reading mesh, order of declaration matter.
    read_by_extension(inputname, volume);
    read_by_extension(inputname, surf);
    FacetAttribute<int> flag(surf);

    if (volume.ncells() != 0) {
        Trace::drop_volume(volume, "volume", {});
        CellFacetAttribute<int> cflag(volume, -1);
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
        transfer_to_surf(volume, cflag, surf, flag);
    }
    else {
        if (surf.nfacets() == 0) {
            std::cerr << "Mesh given: " << inputname << std::endl << "is empty, we cannot compute deformation." << std::endl;
            return 1;
        }
        Trace::drop_surface(surf, "surface", {});
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


    // RUNNING
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
    //disp_polycube(m, U, "cubecover");

    correct_param(m, cfflags, U);

    remove_outerbox(m, cfflags, inside);
    disp_polycube(m, U, cfflags, "corrected");


    // SAVING

    write_by_extension(remeshedname, m, { {{"U", U.ptr}},{},{{"flag", cfflags.ptr}},{} });

    std::ofstream ofs(remeshedflaggingname);
    if (!ofs.is_open()) {
        std::cerr << "Failed opening of flags at : " << remeshedflaggingname << std::endl;
        abort();
    }
    FOR(cf, m.ncells() * 4) ofs << cfflags[cf] << '\n';
    ofs.close();

    FOR(v, m.nverts()) std::swap(m.points[v], U[v]);
    write_by_extension(polycuboidname, m, { {},{},{{"flag", cfflags.ptr}},{} });


    Trace::conclude();
    return 0;
}