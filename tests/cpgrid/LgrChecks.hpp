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


#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>


namespace Opm
{

void checkReferenceElemParentCellVolume(Dune::cpgrid::HierarchicIterator it,
                                        const Dune::cpgrid::HierarchicIterator& endIt);

void checkReferenceElemParentCellCenter(Dune::cpgrid::HierarchicIterator it,
                                        const Dune::cpgrid::HierarchicIterator& endIt,
                                        double expected_total_children);

void checkGlobalCellBounds(const Dune::CpGrid& grid,
                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data,
                           bool lgrsHaveBlockShape,
                           bool isGlobalRefined);

void checkGlobalCellBounds(const std::vector<int>& globalCell,
                           const std::array<int, 3>& logicalCartesianSize);

void checkGlobalCellBoundsConsistencyLevelZeroAndLeaf(const std::vector<int>& globalCell_l0,
                                                      const std::vector<int>& globalCell_leaf);

void checkFatherAndSiblings(const Dune::cpgrid::Entity<0>& element,
                            int preAdaptMaxLevel,
                            double expected_total_children,
                            const Dune::CpGrid& grid);

void checkGridBasicHiearchyInfo(const Dune::CpGrid& grid,
                                const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                int preAdaptMaxLevel = 0);

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

void checkVertexGlobalIds(const Dune::CpGrid& grid,
                          int expected_vertex_ids,
                          int levelOrLeaf);

void checkConsecutiveChildGlobalIdsPerParent(const Dune::CpGrid& grid);

void checkEqualVertexFaceCellSizes(const Dune::CpGrid& grid,
                                   const Dune::CpGrid& other_grid);

void checkEqualVertexGeometrySet(const Dune::CpGrid& grid,
                                 const Dune::CpGrid& other_grid);

void checkEqualCellGeometrySet(const Dune::CpGrid& grid,
                               const Dune::CpGrid& other_grid);

template<typename GridView>
void checkEqualIntersectionsGeometry(const GridView& grid_view,
                                     const Dune::cpgrid::Entity<0>& element,
                                     const GridView& equiv_grid_view,
                                     const Dune::cpgrid::Entity<0>& other_grid_elem);

void checkLeafGridGeometryEquality(const Dune::CpGrid& grid,
                                   const Dune::CpGrid& other_grid);

template<typename T>
bool areClose(const T& cont1,
              const T& cont2);

void adaptGridWithParams(Dune::CpGrid& grid,
                         const std::array<int,3>& cells_per_dim,
                         const std::vector<int>& markedCells);

void adaptGrid(Dune::CpGrid& grid,
               const std::vector<int>& markedCells);

void checkGridWithLgrs(const Dune::CpGrid& grid,
                       const std::vector<std::array<int,3>>& cells_per_dim_vec,
                       const std::vector<std::string>& lgr_name_vec,
                       bool gridHasBeenGlobalRefined = false);

} // namespace Opm

void Opm::checkReferenceElemParentCellVolume(Dune::cpgrid::HierarchicIterator it,
                                             const Dune::cpgrid::HierarchicIterator& endIt)
{
    double reference_elem_parent_cell_volume = std::accumulate(it, endIt, 0.0,
                                                               [](double sum, const Dune::cpgrid::Entity<0>& child) {
                                                                   return sum + child.geometryInFather().volume();
                                                               });
    BOOST_CHECK_CLOSE(reference_elem_parent_cell_volume, 1.0, 1e-13);
}

void Opm::checkReferenceElemParentCellCenter(Dune::cpgrid::HierarchicIterator it,
                                             const Dune::cpgrid::HierarchicIterator& endIt,
                                             double expected_total_children)
{
    int total_children = 0;
    std::set<int> levels{};
    std::array<double,3> reference_elem_parent_cell_center = {0.0, 0.0, 0.0}; // Expected {.5,.5,.5}

    for (; it != endIt; ++it) {
        levels.insert(it->level());
        for (int c = 0; c < 3; ++c) {
            reference_elem_parent_cell_center[c] += it-> geometryInFather().center()[c];
        }
        ++total_children;
    }
    BOOST_CHECK_EQUAL(total_children, expected_total_children);
    BOOST_CHECK_EQUAL( levels.size(), 1); // all children belong to the same level.

    for (int c = 0; c < 3; ++c) {
        reference_elem_parent_cell_center[c]/= total_children;
        BOOST_CHECK_CLOSE(reference_elem_parent_cell_center[c], .5, 1e-13);
    }
}

void Opm::checkGlobalCellBounds(const Dune::CpGrid& grid,
                                const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data,
                                bool lgrsHaveBlockShape,
                                bool isGlobalRefined)
{
    const int maxLevel = grid.maxLevel();
    for (int level = 0; level < maxLevel; ++level) { // level == maxLevel -> leaf grid view
        const auto& level_data = data[level];

        if (level_data->globalCell().empty()) // If level>0, in some processes, LGR might be empty.
            continue;
        checkGlobalCellBounds(level_data->globalCell(), level_data->logicalCartesianSize());
    }

    if(lgrsHaveBlockShape && !isGlobalRefined) {
        checkGlobalCellBoundsConsistencyLevelZeroAndLeaf(data.front()->globalCell(), data.back()->globalCell());
    }
    else {
        checkGlobalCellBounds(data.back()->globalCell(), data.back()->logicalCartesianSize());
    }
}

void Opm::checkGlobalCellBounds(const std::vector<int>& globalCell,
                                const std::array<int, 3>& logicalCartesianSize)
{
    const auto [itMin, itMax] = std::minmax_element(globalCell.begin(), globalCell.end());
    BOOST_CHECK( *itMin >= 0);
    // An LGR can have cells distributed across different processes, so the minimum cell global id may not be zero in each process.

    const auto maxCartesianIdxLevel = logicalCartesianSize[0]*logicalCartesianSize[1]*logicalCartesianSize[2];
    BOOST_CHECK( *itMax < maxCartesianIdxLevel);
}

void Opm::checkGlobalCellBoundsConsistencyLevelZeroAndLeaf(const std::vector<int>& globalCell_l0,
                                                           const std::vector<int>& globalCell_leaf)
{
    const auto [itMinL0, itMaxL0] = std::minmax_element(globalCell_l0.begin(), globalCell_l0.end());
    const auto [itMinLeaf, itMaxLeaf] = std::minmax_element(globalCell_leaf.begin(), globalCell_leaf.end());
    BOOST_CHECK_EQUAL( *itMinL0, *itMinLeaf);
    BOOST_CHECK_EQUAL( *itMaxL0, *itMaxLeaf);
}

void Opm::checkFatherAndSiblings(const Dune::cpgrid::Entity<0>& element,
                                 int preAdaptMaxLevel,
                                 double expected_total_children,
                                 const Dune::CpGrid& grid)
{
    const auto& father = element.father();
    const auto& fatherLevelData = grid.currentData()[father.level()];

    BOOST_CHECK( father.level() <= preAdaptMaxLevel );
    BOOST_CHECK_EQUAL( fatherLevelData ->getMark(father), 1);
    BOOST_CHECK( !father.isLeaf()); // Father vanished during refinement.
    BOOST_CHECK( father.mightVanish() );
    BOOST_CHECK( father == element.getOrigin() );
    /** TODO: Refactor getOrigin() to always return the level-zero ancestor,
        even when the grid has undergone multiple refinements. */

    auto itFather = father.hbegin(grid.maxLevel());
    const auto& endItFather = father.hend(grid.maxLevel());
    // If itFather != endItFather and !father.isLeaf() (if dristibuted_data_ is empty).
    BOOST_CHECK( itFather != endItFather );
    checkReferenceElemParentCellVolume( itFather, endItFather);
    checkReferenceElemParentCellCenter(itFather, endItFather, expected_total_children);

    const auto& [child_level, siblings_list] = fatherLevelData->getChildrenLevelAndIndexList(father.index());

    BOOST_CHECK_EQUAL( child_level, element.level());
    BOOST_CHECK(std::any_of(siblings_list.begin(), siblings_list.end(),
                            [&](int sibling) { return sibling == element.index(); }));
    BOOST_CHECK_EQUAL( siblings_list.size(), expected_total_children);
}

void Opm::checkGridBasicHiearchyInfo(const Dune::CpGrid& grid,
                                     const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                     int preAdaptMaxLevel)
{
    const int maxLevel = grid.maxLevel(); // Leaf Grid View has index maxLevel +1

    for (int level = preAdaptMaxLevel+1; level <= maxLevel+1; ++level) {
        bool isLeaf = (level == maxLevel +1);
        const auto& elements = isLeaf? Dune::elements(grid.leafGridView()) : Dune::elements(grid.levelGridView(level));
        for (const auto& element : elements) {
            BOOST_CHECK_EQUAL( element.hasFather(), static_cast<bool>(element.level())); // level>0 -> 1 (true); level=0 (false)
            BOOST_CHECK( element.getOrigin().level() <= preAdaptMaxLevel);

            if (isLeaf) {
                BOOST_CHECK( !element.mightVanish() );
                BOOST_CHECK( grid.getMark(element) == 0); // postAdapt() has been called, therefore every element gets marked with 0
            }

            auto it = element.hbegin(maxLevel);
            const auto endIt = element.hend(maxLevel);

            if (element.hasFather() && (element.level()>preAdaptMaxLevel)) {
                // second bool discards leaf element born in refined level grids that had not been marked for refinement in the last adapt call.
                int subdivisionsIdx = element.level()-1-preAdaptMaxLevel;
                const auto expected_total_children = cells_per_dim_vec[subdivisionsIdx][0]*cells_per_dim_vec[subdivisionsIdx][1]*cells_per_dim_vec[subdivisionsIdx][2];
                BOOST_CHECK_CLOSE(element.geometryInFather().volume(), 1./expected_total_children, 1e-6);
                checkFatherAndSiblings(element.getLevelElem(), preAdaptMaxLevel, expected_total_children, grid);

                if (!preAdaptMaxLevel) {
                    // If there is no nested refinement, entity.isLeaf() and it == endIt (if dristibuted_data_ is empty).
                    BOOST_CHECK( it == endIt);
                    BOOST_CHECK( element.isLeaf() );
                    BOOST_CHECK( element.isNew() );
                    BOOST_CHECK( !element.mightVanish() );
                }
            }
            else {
                if (preAdaptMaxLevel) {
                    BOOST_CHECK( element.level() <= preAdaptMaxLevel);
                }
                else {
                    BOOST_CHECK_EQUAL( element.level(), 0); // If there is no nested refinement.
                    BOOST_CHECK( element.getOrigin() == element.getLevelElem() );
                    BOOST_CHECK_THROW(element.father(), std::logic_error);
                    BOOST_CHECK_THROW(element.geometryInFather(), std::logic_error);
                    BOOST_CHECK( !element.isNew() );
                }
            }
        }
    }

    /** TODO: Update isNew tag when refine. */
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

    if (maxLevel>=1) { // maxLevel == 0 -> level zero and leaf grids are equal
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
            const auto [all_levelOrLeaf_cell_ids, displLevelOrLeaf] = allGatherv(levelOrLeaf_interior_cell_global_ids, grid.comm());

            const std::unordered_set<int> all_levelOrLeaf_cell_ids_set(all_levelOrLeaf_cell_ids.begin(), all_levelOrLeaf_cell_ids.end());

            // Assert the expected size and uniqueness.
            BOOST_CHECK_EQUAL(static_cast<int>(all_levelOrLeaf_cell_ids.size()), levelOrLeaf_global_cell_count);
            BOOST_CHECK_EQUAL(all_levelOrLeaf_cell_ids.size(), all_levelOrLeaf_cell_ids_set.size());
        }
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

    if (maxLevel>=1) { // maxLevel == 0 -> level zero and leaf grids are equal
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
}

void Opm::checkLgrNameToLevel(const Dune::CpGrid& grid,
                              const std::vector<std::string>& lgr_name_vec)
{
    if (grid.maxLevel()) {
        const auto& lgrNameToLevel =  grid.getLgrNameToLevel();
        for (const auto& [lgrName, level] : lgrNameToLevel) {
            if (level){ // If the grid has been refined only once.
                BOOST_CHECK_EQUAL( lgrName, lgr_name_vec[level-1]);
            }
            else {
                BOOST_CHECK_EQUAL( lgrName, "GLOBAL");
            }
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
    Parser parser;
    const auto deck = parser.parseString(deck_string);
    EclipseState ecl_state(deck);
    EclipseGrid eclipse_grid = ecl_state.getInputGrid();

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
    auto [allGlobalIds_verts, displVertex ] = allGatherv(localVertexIds_vec, grid.comm());
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
template<typename T>
bool Opm::areClose(const T& cont1, const T& cont2)
{
    return (std::abs(cont1[0] - cont2[0]) < 1e-12) && (std::abs(cont1[1] - cont2[1]) < 1e-12) && (std::abs(cont1[2] - cont2[2])< 1e-12);
}

void Opm::checkEqualVertexFaceCellSizes(const Dune::CpGrid& grid, const Dune::CpGrid& other_grid)
{
    BOOST_CHECK_EQUAL(grid.size(3), other_grid.size(3));
    BOOST_CHECK_EQUAL(grid.numFaces(), other_grid.numFaces());
    BOOST_CHECK_EQUAL(grid.size(0), other_grid.size(0));
}

void Opm::checkEqualVertexGeometrySet(const Dune::CpGrid& grid, const Dune::CpGrid& other_grid)
{
    const auto& grid_vertices =  Dune::vertices(grid.leafGridView());
    const auto& other_grid_vertices =  Dune::vertices(other_grid.leafGridView());

    for(const auto& grid_vertex : grid_vertices) {
        // find matching vertex (needed as ordering is allowed to be different
        bool matching_vertex_found = false;
        for (const auto& other_grid_vertex : other_grid_vertices) {
            if (!areClose(grid_vertex.geometry().center(), other_grid_vertex.geometry().center() ))
                continue;
            for(const auto& coord: grid_vertex.geometry().center()) {
                BOOST_TEST(std::isfinite(coord));
            }
            matching_vertex_found = true;
        }
        BOOST_CHECK(matching_vertex_found);
    }
}

void Opm::checkEqualCellGeometrySet(const Dune::CpGrid& grid, const Dune::CpGrid& other_grid)
{
    const auto& grid_view = grid.leafGridView();
    const auto& equiv_grid_view = other_grid.leafGridView();

    const auto& grid_elements =  Dune::elements(grid_view);
    const auto& other_grid_elements =  Dune::elements(equiv_grid_view);

    for(const auto& element : grid_elements) {
        // find matching element (needed as ordering is allowed to be different
        bool matching_elem_found = false;
        for (const auto& other_grid_elem : other_grid_elements) {
            if (!areClose(element.geometry().center(), other_grid_elem.geometry().center() ))
                continue;
            matching_elem_found = true;

            const auto& elem_geo = element.geometry();
            for(const auto& coord : elem_geo.center()) {
                BOOST_TEST(std::isfinite(coord));
            }

            BOOST_CHECK_CLOSE(elem_geo.volume(), other_grid_elem.geometry().volume(), 1e-8);

            const int elemNumFaces = grid.numCellFaces(element.index());
            const int equivElemNumFaces = other_grid.numCellFaces(other_grid_elem.index());
            BOOST_CHECK_EQUAL( elemNumFaces, equivElemNumFaces);

            // For coarse cells with more than 6 faces, this is the case when they touch lgr boundaries, finding matching
            // intersections is not implemented yet.
            if (elemNumFaces>6) {
                continue;
            }
            checkEqualIntersectionsGeometry(grid_view, element, equiv_grid_view, other_grid_elem);
            BOOST_CHECK(matching_elem_found);
        }
    }
}

template<typename GridView>
void Opm::checkEqualIntersectionsGeometry(const GridView& grid_view,
                                          const Dune::cpgrid::Entity<0>& element,
                                          const GridView& equiv_grid_view,
                                          const Dune::cpgrid::Entity<0>& other_grid_elem)
{
    const auto& elemIntersections = Dune::intersections(grid_view, element);
    for(const auto& intersection: elemIntersections) {
        // find matching intersection (needed as ordering is allowed to be different
        bool matching_intersection_found = false;
        const auto& equivElemIntersections = Dune::intersections(equiv_grid_view, other_grid_elem);
        for(const auto& intersection_match : equivElemIntersections) {
            if(intersection_match.indexInInside() == intersection.indexInInside()) {
                BOOST_CHECK(intersection_match.neighbor() == intersection.neighbor());

                if(intersection.neighbor()) {
                    BOOST_CHECK(intersection_match.indexInOutside() == intersection.indexInOutside());
                }

                BOOST_CHECK( areClose(intersection_match.centerUnitOuterNormal(), intersection.centerUnitOuterNormal()) );

                const auto& geom_match = intersection_match.geometry();
                const auto& geom =  intersection.geometry();
                bool closeGeomCenter = areClose(geom_match.center(), geom.center());
                if (!closeGeomCenter) {
                    break; // Check next intersection_match
                }

                BOOST_CHECK_CLOSE(geom_match.volume(), geom.volume(), 1e-6);
                BOOST_CHECK(areClose(geom_match.center(), geom.center()));
                BOOST_CHECK(geom_match.corners() == geom.corners());

                decltype(geom.corner(0)) sum_match{}, sum{};

                for(int cor = 0; cor < geom.corners(); ++cor) {
                    sum += geom.corner(cor);
                    sum_match += geom_match.corner(1);
                }
                BOOST_CHECK(areClose(sum, sum_match));
                matching_intersection_found = true;
                break;
            }
        }
        BOOST_CHECK(matching_intersection_found);
    }
}

void Opm::checkLeafGridGeometryEquality(const Dune::CpGrid& grid, const Dune::CpGrid& other_grid)
{
    checkEqualVertexFaceCellSizes(grid, other_grid);
    checkEqualVertexGeometrySet(grid, other_grid);
    checkEqualCellGeometrySet(grid, other_grid);
}

void Opm::adaptGridWithParams(Dune::CpGrid& grid,
                              const std::array<int,3>& cells_per_dim,
                              const std::vector<int>& markedCells)
{
    const int startingGridIdx = grid.currentData().size() -1; // size before calling adapt
    std::vector<int> assignRefinedLevel(grid.currentData()[startingGridIdx]->size(0));

    for (const auto& elemIdx : markedCells)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(grid.currentData()[startingGridIdx]), elemIdx, true);
        grid.mark(1, elem);
        assignRefinedLevel[elemIdx] = grid.maxLevel() + 1;
        BOOST_CHECK( grid.getMark(elem) == 1);
        BOOST_CHECK( elem.mightVanish() == true);
    }
    grid.preAdapt();
    grid.adapt({cells_per_dim}, assignRefinedLevel, {"LGR"+std::to_string(grid.maxLevel() +1)});
    grid.postAdapt();
}

void Opm::adaptGrid(Dune::CpGrid& grid,
                    const std::vector<int>& markedCells)
{
    const auto& leafGridView = grid.currentData().back();
    for (const auto& elemIdx : markedCells)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*leafGridView, elemIdx, true);
        grid.mark(1, elem);
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();
}

void Opm::checkGridWithLgrs(const Dune::CpGrid& grid,
                            const std::vector<std::array<int,3>>& cells_per_dim_vec,
                            const std::vector<std::string>& lgr_name_vec,
                            bool gridHasBeenGlobalRefined)
{
    const auto& data = grid.currentData();

    checkLgrNameToLevel(grid, lgr_name_vec);
    checkVertexAndFaceIndexAreNonNegative(grid);
    checkFaceHas4VerticesAndMax2NeighboringCells(grid, data);
    checkLocalIndicesMatchMapper(grid); // Decide if it's worth to keep it
    checkGridBasicHiearchyInfo(grid, cells_per_dim_vec);

    checkCellGlobalIdUniquenessForInteriorCells(grid, data);
    /** Vertex global id uniqueness can't be guaranteed because the overlap layer size is set to 1.
        This means cells sharing only corners or edges (not faces) with interior cells aren't visible to the process,
        potentially causing "multiple ids" for the same vertices.
        To achieve unique vertex ids in some cases, use loadBalance(parts, false, true),
        where 'parts' sets cell ranks, 'false' means ownerFirst, and 'true' adds corner cells.
        The overlapLayerSize defaults to 1.
    */
    checkGridLocalAndGlobalIdConsistency(grid, data);

    checkGlobalCellBounds(grid, data, /* lgrsHaveBlockShape */ true, gridHasBeenGlobalRefined);
}

#endif // OPM_LGRCHECKS_HEADER_INCLUDED

