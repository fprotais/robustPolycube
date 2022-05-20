#include <flagging/naive_flagging.h>
#include <utils/trace.h>
#include <filesystem>

#ifndef DEBUG_GRAPHITE_PATH
#define DEBUG_GRAPHITE_PATH "C:/fprotais/softwares/graphite/build/Windows/bin/Release/graphite.exe"
#endif 

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)


int main(int argc, char** argv) {
	Trace::initialize(DEBUG_GRAPHITE_PATH);
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " mesh.ext flaggingfile" << std::endl;
        std::cerr << "Input:" << std::endl;
        std::cerr << "- mesh.ext must contain a tet mesh or tri mesh" << std::endl;
        std::cerr << "Output: [OPTIONAL]" << std::endl;
        std::cerr << "- flaggingfile is were we save flagging for mesh, see README.md for details, default name will be 'mesh.flags'" << std::endl;
        std::cerr << std::endl;
        std::cerr << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cerr << "For more details, see README." << std::endl;
        std::cerr << "contact: francoisprotais@gmail.com" << std::endl;
        return 1;
    }
    std::string filename(argv[1]);
    std::string savename = std::filesystem::path(filename).stem().string() + ".flags";
    if (argc > 2) savename = argv[2];

    Tetrahedra m;
    read_by_extension(filename, m);

    if (m.ncells() != 0) {
        Trace::drop_volume(m, "volume", {});
        CellFacetAttribute<int> flag(m, -1);
        generate_naive_flagging(m, flag);
        Trace::drop_cellfacet_scalar(m, flag, "flags", -1, true);
        std::ofstream ofs(savename);
        if (!ofs.is_open()) {
            std::cerr << "Failed opening of flags at : " << savename << std::endl;
            abort();
        }
        FOR(cf, m.ncells() * 4) ofs << flag[cf] << '\n';
        ofs.close();
    }
    else {
        Triangles surf;
        read_by_extension(filename, surf);
        if (surf.nfacets() == 0) {
            std::cerr << "Mesh given: " << filename << std::endl << "is empty, we cannot compute flagging." << std::endl;
            return 1;
        }
        Trace::drop_surface(surf, "surface", {});
        FacetAttribute<int> flag(surf);
        generate_naive_flagging(surf, flag);
        Trace::drop_facet_scalar(surf, flag, "flagging");
        std::ofstream ofs(savename);
        if (!ofs.is_open()) {
            std::cerr << "Failed opening of flags at : " << savename << std::endl;
            abort();
        }
        FOR(f, surf.nfacets()) ofs << flag[f] << '\n';
        ofs.close();
    }

	Trace::conclude();
    return 0;
}