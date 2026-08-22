#include <config.h>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

#if HAVE_OPM_COMMON
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#endif

#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

#if HAVE_OPM_COMMON
Opm::ParseContext makeFlowLikeParseContext()
{
    return Opm::ParseContext({
        {Opm::ParseContext::PARSE_RANDOM_SLASH, Opm::InputErrorAction::IGNORE},
        {Opm::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Opm::InputErrorAction::WARN},
        {Opm::ParseContext::SUMMARY_UNKNOWN_WELL, Opm::InputErrorAction::WARN},
        {Opm::ParseContext::SUMMARY_UNKNOWN_GROUP, Opm::InputErrorAction::WARN},
    });
}
#endif

void writeGridFromCpGrid(const Dune::CpGrid& grid, const std::string& output_file)
{
    const int ndims = 3;
    const int ncells = grid.numCells();
    const int nfaces = grid.numFaces();
    const int nnodes = grid.numVertices();
    const int ncellfaces = grid.numCellFaces();

    int nfacenodes = 0;
    std::vector<int> face_nodepos;
    std::vector<int> face_nodes;
    face_nodepos.reserve(nfaces + 1);
    face_nodes.reserve(4 * nfaces);
    face_nodepos.push_back(0);
    for (int f = 0; f < nfaces; ++f) {
        const int nv = grid.numFaceVertices(f);
        for (int lv = 0; lv < nv; ++lv) {
            face_nodes.push_back(grid.faceVertex(f, lv));
        }
        nfacenodes += nv;
        face_nodepos.push_back(nfacenodes);
    }

    std::vector<int> face_cells;
    face_cells.reserve(2 * nfaces);
    for (int f = 0; f < nfaces; ++f) {
        face_cells.push_back(grid.faceCell(f, 0));
        face_cells.push_back(grid.faceCell(f, 1));
    }

    std::vector<int> cell_facepos;
    std::vector<int> cell_faces;
    std::vector<int> cell_facetag;
    cell_facepos.reserve(ncells + 1);
    cell_faces.reserve(ncellfaces);
    cell_facetag.reserve(ncellfaces);
    cell_facepos.push_back(0);
    const auto c2f = Opm::UgGridHelpers::cell2Faces(grid);
    for (int c = 0; c < ncells; ++c) {
        const auto row = c2f[c];
        for (auto it = row.begin(); it != row.end(); ++it) {
            const int face = *it;
            cell_faces.push_back(face);

            // Export the exact CpGrid face tag used by transmissibility logic
            const int tag = Opm::UgGridHelpers::faceTag(grid, it);
            cell_facetag.push_back(tag);
        }
        cell_facepos.push_back(static_cast<int>(cell_faces.size()));
    }

    if (static_cast<int>(cell_faces.size()) != ncellfaces) {
        throw std::runtime_error("Internal inconsistency: numCellFaces() total mismatch");
    }

    std::ofstream out(output_file);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + output_file);
    }

    out << std::setprecision(17) << std::fixed;

    const int has_tag = 1;
    const int has_indexmap = 1;
    out << ndims << ' ' << ncells << ' ' << nfaces << ' '
        << nnodes << ' ' << nfacenodes << ' ' << ncellfaces << ' '
        << has_tag << ' ' << has_indexmap << '\n';

    const auto cartdims = grid.logicalCartesianSize();
    out << cartdims[0] << ' ' << cartdims[1] << ' ' << cartdims[2] << '\n';

    for (int v = 0; v < nnodes; ++v) {
        const auto& p = grid.vertexPosition(v);
        out << p[0] << ' ' << p[1] << ' ' << p[2] << ' ';
    }
    out << '\n';

    for (const int value : face_nodepos) {
        out << value << ' ';
    }
    out << '\n';

    for (const int value : face_nodes) {
        out << value << ' ';
    }
    out << '\n';

    for (const int value : face_cells) {
        out << value << ' ';
    }
    out << '\n';

    for (int f = 0; f < nfaces; ++f) {
        out << grid.faceArea(f) << ' ';
    }
    out << '\n';

    for (int f = 0; f < nfaces; ++f) {
        const auto& c = grid.faceCentroid(f);
        out << c[0] << ' ' << c[1] << ' ' << c[2] << ' ';
    }
    out << '\n';

    for (int f = 0; f < nfaces; ++f) {
        const auto& n = grid.faceNormal(f);
        out << n[0] << ' ' << n[1] << ' ' << n[2] << ' ';
    }
    out << '\n';

    for (const int value : cell_facepos) {
        out << value << ' ';
    }
    out << '\n';

    for (std::size_t i = 0; i < cell_faces.size(); ++i) {
        out << cell_faces[i] << ' ' << cell_facetag[i] << ' ';
    }
    out << '\n';

    const auto& global_cell = grid.globalCell();
    for (const int value : global_cell) {
        out << value << ' ';
    }
    out << '\n';

    for (int c = 0; c < ncells; ++c) {
        out << grid.cellVolume(c) << ' ';
    }
    out << '\n';

    for (int c = 0; c < ncells; ++c) {
        const auto& cc = grid.cellCentroid(c);
        out << cc[0] << ' ' << cc[1] << ' ' << cc[2] << ' ';
    }
    out << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

#if !HAVE_OPM_COMMON
    std::cerr << "export_grid requires HAVE_OPM_COMMON to parse .DATA files." << std::endl;
    return 2;
#else
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <case.DATA> <output.grid>" << std::endl;
        return 1;
    }

    const std::string data_file = argv[1];
    const std::string output_file = argv[2];

    try {
        Dune::CpGrid cpgrid;
        bool built_from_data = false;

        try {
            Opm::Parser parser;
            const auto deck = parser.parseFile(data_file, makeFlowLikeParseContext());

            Opm::EclipseGrid eclipse_grid(deck);
            Opm::EclipseState ecl_state(deck);

            cpgrid.processEclipseFormat(&eclipse_grid, &ecl_state, false);
            built_from_data = true;
        } catch (const std::exception& e) {
            const std::filesystem::path data_path(data_file);
            const std::filesystem::path egrid_path =
                data_path.parent_path() / (data_path.stem().string() + ".EGRID");
            if (!std::filesystem::exists(egrid_path)) {
                throw std::runtime_error(
                    std::string("Failed to parse DATA and no EGRID fallback found. DATA error: ") + e.what());
            }

            std::cerr << "Warning: failed to build grid directly from DATA (" << e.what() << ")\n"
                      << "Falling back to EGRID: " << egrid_path << "\n";
            Opm::EclipseGrid eclipse_grid(egrid_path.string());
            cpgrid.processEclipseFormat(&eclipse_grid, nullptr, false);
        }

        writeGridFromCpGrid(cpgrid, output_file);

        std::cout << "Wrote grid: " << output_file << "\n"
                  << "  source             : " << (built_from_data ? "DATA" : "EGRID fallback") << "\n"
                  << "  cells/faces/vertices/cellFaces = "
                  << cpgrid.numCells() << " / "
                  << cpgrid.numFaces() << " / "
                  << cpgrid.numVertices() << " / "
                  << cpgrid.numCellFaces() << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Failed to export grid from DATA: " << e.what() << std::endl;
        return 3;
    }

    return 0;
#endif
}
