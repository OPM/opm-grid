/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_LGRCHECKS_HEADER_INCLUDED
#define OPM_LGRCHECKS_HEADER_INCLUDED

#include <boost/test/unit_test.hpp>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/common/version.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>


#include <array>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>



namespace Opm
{

void checkAverageCenterAndVolume(Dune::cpgrid::HierarchicIterator it,
                                 const Dune::cpgrid::HierarchicIterator& endIt,
                                 double total_children);

void checkGlobalCellBounds(const Dune::CpGrid& grid,
                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data);

void checkGlobalCellBounds(const std::vector<int>& globalCell,
                           const std::array<int, 3>& logicalCartesianSize);

void checkEqMinMaxGlobalCellLevelZeroAndLeaf(const std::vector<int>& globalCell_l0,
                                             const std::vector<int>& globalCell_leaf);

void checkFatherAndSiblings(const Dune::cpgrid::Entity<0>& element,
                            int level,
                            double expected_total_children,
                            const Dune::CpGrid& grid);

void checkGridBasicHiearchyInfo(const Dune::CpGrid& grid,
                                const std::vector<std::array<int,3>>& cells_per_dim_vec);

void checkLocalIndicesMatchMapper(const Dune::CpGrid& grid);

void checkCellGlobalIdUniquenessForInteriorCells(const Dune::CpGrid& grid,
                                                 const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data);

void checkGridLocalAndGlobalIdConsistency(const Dune::CpGrid& grid,
                                          const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data);

void checkLgrNameToLevel(const Dune::CpGrid& grid,
                         const std::vector<std::string>& lgr_name_vec);

void checkVertexAndFaceIndexAreNonNegative(const Dune::CpGrid& grid);

void checkFaceHas4VerticesAndMax2NeighboringCells(const Dune::CpGrid& grid,
                                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data);

void createGridAndAddLgrs(Dune::CpGrid& grid,
                          const std::string& deck_string,
                          const std::vector<std::array<int, 3>>& cells_per_dim_vec,
                          const std::vector<std::array<int, 3>>& startIJK_vec,
                          const std::vector<std::array<int, 3>>& endIJK_vec,
                          const std::vector<std::string>& lgr_name_vec);

void createGridAndAddLgrs(Dune::CpGrid& grid,
                          const std::array<double, 3>& cell_sizes,
                          const std::array<int, 3>& grid_dim,
                          const std::vector<std::array<int, 3>>& cells_per_dim_vec,
                          const std::vector<std::array<int, 3>>& startIJK_vec,
                          const std::vector<std::array<int, 3>>& endIJK_vec,
                          const std::vector<std::string>& lgr_name_vec);

void checkExpectedVertexGlobalIdsCount(const Dune::CpGrid& grid,
                                       const std::vector<int>& expected_vertex_ids_per_lgr,
                                       int leaf_expected_vertex_ids);

void checkVertexGlobalIds(const Dune::CpGrid& grid, int expected_vertex_ids, int levelOrLeaf);

void checkConsecutiveChildGlobalIdsPerParent(const Dune::CpGrid& grid);


} // namespace Opm




void Opm::checkAverageCenterAndVolume(Dune::cpgrid::HierarchicIterator it,
                                      const Dune::cpgrid::HierarchicIterator& endIt,
                                      double total_children)
{
    double referenceElemOneParent_volume_it{};
    std::array<double,3> referenceElem_entity_center_it; // Expected {.5,.5,.5}

    for (; it != endIt; ++it)
    {
        referenceElemOneParent_volume_it += it-> geometryInFather().volume();
        for (int c = 0; c < 3; ++c)
        {
            referenceElem_entity_center_it[c] += (it-> geometryInFather().center())[c];
        }
    }
    BOOST_CHECK_CLOSE(referenceElemOneParent_volume_it, 1, 1e-13);

    for (int c = 0; c < 3; ++c)
    {
        referenceElem_entity_center_it[c]/= total_children;
        BOOST_CHECK_CLOSE(referenceElem_entity_center_it[c], .5, 1e-13);
    }
}

void Opm::checkGlobalCellBounds(const Dune::CpGrid& grid,
                                const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data)
{
    const int maxLevel = grid.maxLevel();
    for (int level = 0; level <= maxLevel; ++level) {
        if (!(data[level] -> globalCell().empty())) { // If level>0, in some processes, LGR might be empty.
            Opm::checkGlobalCellBounds(data[level]->globalCell(), data[level]->logicalCartesianSize());
        }
    }
    Opm::checkGlobalCellBounds(data.back()->globalCell(), data.back()->logicalCartesianSize());
    Opm::checkEqMinMaxGlobalCellLevelZeroAndLeaf(data.front()->globalCell(), data.back()->globalCell());
}

void Opm::checkGlobalCellBounds(const std::vector<int>& globalCell,
                                const std::array<int, 3>& logicalCartesianSize)
{
    const auto [itMin, itMax] = std::minmax_element(globalCell.begin(), globalCell.end());
    BOOST_CHECK( *itMin >= 0);
    //Note:  An LGR can have cells distributed across different processes, so the minimum cell global id may not be zero in all processes.

    const auto maxCartesianIdxLevel = logicalCartesianSize[0]*logicalCartesianSize[1]*logicalCartesianSize[2];
    BOOST_CHECK( *itMax < maxCartesianIdxLevel);
}

void Opm::checkEqMinMaxGlobalCellLevelZeroAndLeaf(const std::vector<int>& globalCell_l0,
                                                  const std::vector<int>& globalCell_leaf)
{
    const auto [itMinL0, itMaxL0] = std::minmax_element(globalCell_l0.begin(), globalCell_l0.end());
    const auto [itMinLeaf, itMaxLeaf] = std::minmax_element(globalCell_leaf.begin(), globalCell_leaf.end());
    BOOST_CHECK_EQUAL( *itMinL0, *itMinLeaf);
    BOOST_CHECK_EQUAL( *itMaxL0, *itMaxLeaf);
}

void Opm::checkFatherAndSiblings(const Dune::cpgrid::Entity<0>& element,
                                 int level,
                                 double expected_total_children,
                                 const Dune::CpGrid& grid)
{
    BOOST_CHECK( element.father().level() == 0); // If there is no nested refinement
    BOOST_CHECK( element.father().isLeaf() == false); // Father vanished during refinement.
    BOOST_CHECK( element.father() == element.getOrigin() );

    auto itFather = element.father().hbegin(grid.maxLevel());
    const auto endItFather = element.father().hend(grid.maxLevel());
    // If itFather != endItFather and !element.father().isLeaf() (if dristibuted_data_ is empty).
    BOOST_CHECK( itFather != endItFather );
    Opm::checkAverageCenterAndVolume(itFather, endItFather, expected_total_children);

    const auto& [child_level, siblings_list] =
        grid.currentData()[element.father().level()]->getChildrenLevelAndIndexList(element.father().index());

    BOOST_CHECK_EQUAL( child_level, level);

    BOOST_CHECK( (std::find(siblings_list.begin(), siblings_list.end(), element.index() ) != siblings_list.end()));
    BOOST_CHECK_EQUAL( siblings_list.size(), expected_total_children);
}

void Opm::checkGridBasicHiearchyInfo(const Dune::CpGrid& grid,
                                     const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    const int maxLevel = grid.maxLevel(); // Leaf Grid View has index maxLevel +1
    for (int level = 0; level <= maxLevel+1; ++level) {
        bool isLeaf = (level == maxLevel +1);
        const auto& elements = isLeaf? Dune::elements(grid.leafGridView()) : Dune::elements(grid.levelGridView(level));
        for (const auto& element : elements) {
            BOOST_CHECK_EQUAL( element.hasFather(), static_cast<bool>(element.level())); // level>0 -> 1 (true); level=0 (false)
            BOOST_CHECK( element.getOrigin().level() == 0);

            auto it = element.hbegin(maxLevel);
            const auto endIt = element.hend(maxLevel);

            if (element.hasFather()) {
                // If there is no nested refinement, entity.isLeaf() and it == endIt (if dristibuted_data_ is empty).
                BOOST_CHECK( it == endIt);
                BOOST_CHECK( element.isLeaf() );

                const auto expected_total_children = cells_per_dim_vec[element.level()-1][0]*cells_per_dim_vec[element.level()-1][1]*cells_per_dim_vec[element.level()-1][2];
                BOOST_CHECK_CLOSE(element.geometryInFather().volume(), 1./expected_total_children, 1e-6);

                checkFatherAndSiblings(element.getLevelElem(), element.level(), expected_total_children, grid);
            }
            else {
                BOOST_CHECK( element.getOrigin() == element.getLevelElem() );
                BOOST_CHECK_THROW(element.father(), std::logic_error);
                BOOST_CHECK_THROW(element.geometryInFather(), std::logic_error);
                BOOST_CHECK_EQUAL( element.level(), 0); // If there is no nested refinement.
            }
        }
    }
}

void Opm::checkLocalIndicesMatchMapper(const Dune::CpGrid& grid)
{
    const auto& leafView = grid.leafGridView();
    Dune::MultipleCodimMultipleGeomTypeMapper leafMapper(leafView, Dune::mcmgElementLayout());

    for (const auto& element : elements(leafView)) {
        BOOST_CHECK_EQUAL(element.index(), leafMapper.index(element));
    }

    const int maxLevel = grid.maxLevel();
    for (int level = 0; level <= maxLevel; ++level) {
        const auto& levelView = grid.levelGridView(level);
        Dune::MultipleCodimMultipleGeomTypeMapper levelMapper(levelView, Dune::mcmgElementLayout());

        for (const auto& element : elements(levelView)) {
            BOOST_CHECK_EQUAL(element.index(), levelMapper.index(element));
        }
    }
}

void Opm::checkCellGlobalIdUniquenessForInteriorCells(const Dune::CpGrid& grid,
                                                      const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data)
{
    const auto& leafGridView = grid.leafGridView();
    const int maxLevel = grid.maxLevel();

    for (int level = 1; level <= maxLevel +1; ++level) {
        bool isLeaf = (level == maxLevel+1);
        const auto& levelOrLeafData = data[level];
        const auto& levelOrLeafGlobalIdSet = levelOrLeafData->globalIdSet();
        const auto& elements = isLeaf?
            Dune::elements(leafGridView, Dune::Partitions::interior) :
            Dune::elements(grid.levelGridView(level), Dune::Partitions::interior);

        std::vector<int> levelOrLeaf_interior_cell_global_ids;
        levelOrLeaf_interior_cell_global_ids.reserve(levelOrLeafData->size(0));

        for (const auto& element : elements) {
            levelOrLeaf_interior_cell_global_ids.push_back(levelOrLeafGlobalIdSet.id(element));
        }

        const int levelOrLeaf_local_interior_cell_count = levelOrLeaf_interior_cell_global_ids.size();
        const int levelOrLeaf_global_cell_count = grid.comm().sum(levelOrLeaf_local_interior_cell_count);

        // Gather all global IDs across MPI ranks.
        const auto [all_levelOrLeaf_cell_ids, displLevelOrLeaf] = Opm::allGatherv(levelOrLeaf_interior_cell_global_ids, grid.comm());

        const std::unordered_set<int> all_levelOrLeaf_cell_ids_set(all_levelOrLeaf_cell_ids.begin(), all_levelOrLeaf_cell_ids.end());

        // Assert the expected size and uniqueness.
        BOOST_CHECK_EQUAL(static_cast<int>(all_levelOrLeaf_cell_ids.size()), levelOrLeaf_global_cell_count);
        BOOST_CHECK_EQUAL(all_levelOrLeaf_cell_ids.size(), all_levelOrLeaf_cell_ids_set.size());
    }
}

void Opm::checkGridLocalAndGlobalIdConsistency(const Dune::CpGrid& grid,
                                               const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data)
{
    const auto& leafGridView = grid.leafGridView();
    const auto& leafLocalIdSet = data.back()->localIdSet();

    const bool isSerialRun = (grid.comm().size() == 0);

    std::unordered_map<int, Dune::cpgrid::Entity<0>> leafElementLocalId;
    std::unordered_map<int, Dune::cpgrid::Entity<3>> leafVertexLocalId;

    const auto& leafElements = Dune::elements(leafGridView);

    // Maps to verify if a level element or vertex belongs to the leaf grid.
    for (const auto& leafElem : leafElements) {
        leafElementLocalId[leafLocalIdSet.id(leafElem)] = leafElem;
    }
    for (const auto& leafVertex : vertices(leafGridView)) {
        leafVertexLocalId[leafLocalIdSet.id(leafVertex)] = leafVertex;
    }

    const int maxLevel = grid.maxLevel(); // Leaf grid view in level = maxLevel+1
    for (int level = 0; level <= maxLevel+1; ++level) {
        bool isLeaf = (level == maxLevel+1);
        const auto& elements = isLeaf? leafElements : Dune::elements(grid.levelGridView(level));
        const auto& vertices = isLeaf? Dune::vertices(leafGridView) : Dune::vertices(grid.levelGridView(level));
        const auto& localIdSet = data[level]->localIdSet();
        const auto& globalIdSet = data[level]->globalIdSet();
        const auto& indexSet = data[level]->indexSet();

        std::unordered_set<int> levelOrLeafIds_set;
        levelOrLeafIds_set.reserve(data[level]->size(0) + data[level]->size(3));

        for (const auto& element : elements) {
            const int localId = localIdSet.id(element);
            const int globalId = globalIdSet.id(element);

            if (isSerialRun) {
                BOOST_CHECK_EQUAL(localId, globalId);
            }

            levelOrLeafIds_set.insert(localId);

            if (!isLeaf) {
                if (auto it = leafElementLocalId.find(localId); it != leafElementLocalId.end()) {
                    BOOST_CHECK(it->second.geometry().center() == element.geometry().center());
                    BOOST_CHECK(it->second.geometry().volume() == element.geometry().volume());
                    BOOST_CHECK(it->second.getLevelElem() == element);
                }
                const int idx = indexSet.index(element);
                BOOST_CHECK_EQUAL(idx, element.index());
            }

            BOOST_CHECK_EQUAL(grid.globalIdSet().id(element),
                              data[element.level()]->globalIdSet().id(element.getLevelElem()));
        }

        BOOST_CHECK_EQUAL(static_cast<int>(levelOrLeafIds_set.size()), data[level]->size(0));

        for (const auto& vertex : vertices) {
            const int localId = localIdSet.id(vertex);
            const int globalId = globalIdSet.id(vertex);

            if (isSerialRun) {
                BOOST_CHECK_EQUAL(localId, globalId);
            }

            levelOrLeafIds_set.insert(localId);

            if (!isLeaf) {
                if (auto it = leafVertexLocalId.find(localId); it != leafVertexLocalId.end()) {
                    BOOST_CHECK(it->second.geometry().center() == vertex.geometry().center());
                }
            }
        }

        BOOST_CHECK_EQUAL(static_cast<int>(levelOrLeafIds_set.size()), data[level]->size(0) + data[level]->size(3));
    }
}

void Opm::checkLgrNameToLevel(const Dune::CpGrid& grid,
                              const std::vector<std::string>& lgr_name_vec)
{
    BOOST_CHECK_EQUAL( grid.maxLevel(), static_cast<int>(lgr_name_vec.size())); // If there is no nested refinement

    const auto& lgrNameToLevel =  grid.getLgrNameToLevel();
    for (const auto& [lgrName, level] : lgrNameToLevel) {
        if (level){
            BOOST_CHECK_EQUAL( lgrName, lgr_name_vec[level-1]);
        }
        else {
            BOOST_CHECK_EQUAL( lgrName, "GLOBAL");
        }
    }
}

void Opm::checkVertexAndFaceIndexAreNonNegative(const Dune::CpGrid& grid)
{
    const auto& leafGridView = grid.leafGridView();
    for (const auto& element : elements(leafGridView)) {
        BOOST_CHECK_EQUAL( element.subEntities(3), 8);
        for (int i = 0; i < 8; ++i){
            BOOST_CHECK( element.subEntity<3>(i).index() >= 0); // valid vertex index
        }
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            BOOST_CHECK( intersection.id() >= 0); // valid face index
        }
    }
}

void Opm::checkFaceHas4VerticesAndMax2NeighboringCells(const Dune::CpGrid& grid,
                                                       const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data)
{
    /** Add level grid check */

    const int numFaces = grid.numFaces();
    for (int face = 0; face < numFaces; ++face){
        BOOST_CHECK_EQUAL( grid.numFaceVertices(face), 4);
        for (int i = 0; i < 4; ++i) {
            BOOST_CHECK( grid.faceVertex(face, i) >= 0); // valid vertex index
        }
        BOOST_CHECK( data.back()->faceToCellSize(face) < 3); // Max. 2 neighboring cells
    }
}

void Opm::createGridAndAddLgrs(Dune::CpGrid& grid,
                               const std::string& deck_string,
                               const std::vector<std::array<int, 3>>& cells_per_dim_vec,
                               const std::vector<std::array<int, 3>>& startIJK_vec,
                               const std::vector<std::array<int, 3>>& endIJK_vec,
                               const std::vector<std::string>& lgr_name_vec)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();

    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);

    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

void Opm::createGridAndAddLgrs(Dune::CpGrid& grid,
                               const std::array<double, 3>& cell_sizes,
                               const std::array<int, 3>& grid_dim,
                               const std::vector<std::array<int, 3>>& cells_per_dim_vec,
                               const std::vector<std::array<int, 3>>& startIJK_vec,
                               const std::vector<std::array<int, 3>>& endIJK_vec,
                               const std::vector<std::string>& lgr_name_vec)
{
    grid.createCartesian(grid_dim, cell_sizes);
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

void Opm::checkExpectedVertexGlobalIdsCount(const Dune::CpGrid& grid,
                                            const std::vector<int>& expected_vertex_ids_per_lgr,
                                            int leaf_expected_vertex_ids)
{
    for (std::size_t lgr = 1; lgr < expected_vertex_ids_per_lgr.size(); ++lgr) {
        checkVertexGlobalIds(grid, expected_vertex_ids_per_lgr[lgr-1], lgr);
    }
    checkVertexGlobalIds(grid, leaf_expected_vertex_ids, grid.maxLevel()+1);
}

void Opm::checkVertexGlobalIds(const Dune::CpGrid& grid, int expected_vertex_ids, int levelOrLeaf)
{
    std::vector<int> localVertexIds_vec;
    const auto& levelOrLeafData = grid.currentData()[levelOrLeaf];
    localVertexIds_vec.reserve(levelOrLeafData->size(3));

    const int leafGridIdx = grid.maxLevel()+1;
    bool isLeaf = (levelOrLeaf == leafGridIdx);

    const auto& verts = isLeaf? Dune::vertices(grid.leafGridView()) : Dune::vertices(grid.levelGridView(levelOrLeaf));
    // Notice that all partition type points are pushed back.
    // Selecting only interior points does not bring us to the expected value.
    std::transform(verts.begin(), verts.end(), std::back_inserter(localVertexIds_vec),
                   [&is = levelOrLeafData->globalIdSet()](const auto& vertex)
                   { return is.id(vertex); });
    auto [allGlobalIds_verts, displVertex ] = Opm::allGatherv(localVertexIds_vec, grid.comm());
    const std::set<int> allGlobalIds_verts_set(allGlobalIds_verts.begin(), allGlobalIds_verts.end());

    BOOST_CHECK_EQUAL( allGlobalIds_verts_set.size(), expected_vertex_ids);
}

void Opm::checkConsecutiveChildGlobalIdsPerParent(const Dune::CpGrid& grid)
{
    const auto& data = grid.currentData(); // data_ or distributed_data_
    const auto& globalIdSetBeforeLoadBalance = grid.globalIdSet();
    const auto& elements = Dune::elements(grid.levelGridView(0));
    const auto& parentToChildrenBeforeLoadBalance = data[0]->getParentToChildren();

    for (const auto& element : elements) {
        const auto& [level, children] = parentToChildrenBeforeLoadBalance[element.index()];

        if (!children.empty()) {
            const auto& levelData = *data[level];
            const std::size_t numChildren = children.size();

            const auto& child0_elem = Dune::cpgrid::Entity<0>(levelData, children[0], true);
            const auto child0_globalId = globalIdSetBeforeLoadBalance.id(child0_elem);

            for (std::size_t idx = 1; idx < numChildren; ++idx) {
                const auto& child_elem = Dune::cpgrid::Entity<0>(levelData, children[idx], true);
                const auto child_globalId = globalIdSetBeforeLoadBalance.id(child_elem);

                BOOST_REQUIRE(numChildren > 1);
                BOOST_CHECK_EQUAL(child0_globalId + idx, child_globalId);
            }
        }
    }
}

#endif // OPM_LGRCHECKS_HEADER_INCLUDED

