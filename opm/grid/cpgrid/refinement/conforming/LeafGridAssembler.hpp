/*
  Copyright 2026 Equinor ASA.

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
#ifndef OPM_GRID_REFINEMENT_LEAFGRIDASSEMBLER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_LEAFGRIDASSEMBLER_HEADER_INCLUDED

#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/OrientedEntityTable.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <opm/common/utility/SparseTable.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

namespace Dune
{
namespace cpgrid
{
class CpGridData;
}
}

namespace Opm
{
namespace Refinement
{

/// Merge level zero and the refined level grids into the leaf view.
///
/// Leaf cells are the unrefined level-zero cells plus all refined cells
/// (children replace their parent at its position in the cell ordering).
/// Parent faces on block boundaries are replaced by the refined mosaic,
/// pairing refined cells with their level-zero neighbors. Corners follow
/// one shared pool: the leaf corner vector starts with all
/// level-zero corners (indices preserved), refined corners that coincide
/// with parent-lattice corners are identified with the level-zero corner
/// through the parent cells (so faults inside blocks identify per cell,
/// not per lattice key), and the rest are appended.
///
/// Ids work through the existing delegation machinery: levels_.leaf_to_level_cells
/// and levels_.corner_history are populated so IdSet routes every leaf entity to
/// the id of its birth-level entity.
///
/// Block-boundary parent faces crossed by a fault are rebuilt from the
/// corner-point processor (split faces against the coarse neighbour(s)); the
/// parent description is needed for that.
///
/// A nested box (one refining another box's level grid) replaces the cells
/// it refines by its own: its corners on its parents' corners are theirs, and
/// its sides face the parent LGR's cells. It must lie inside its parent LGR,
/// away from the parent's boundary, one level deep, with no fault or pinch-out
/// on its sides.
///
/// The constructor stores the input and computes what no stage changes: the box
/// tree and the basic level-zero lookups, held as const members. build() runs the
/// stages below, one private method each, in this order. Each stage writes the
/// members listed with it; a later stage reads them through the const accessors.
///  1. collectBoxData: boxes_ (except their cornerToLeaf and faceToLeaf), mask0_
///     and nestedBoxOf_.
///  2. buildCornerPool: leaf_ with its corners, leafCornerOrigin_,
///     numLeafCorners_, and each box's cornerToLeaf.
///  3. buildLeafCells: leafCells_, leafIdxOfCell0_, leafIdxOfLevelCell_,
///     level0ToChildren_, numLeafCells_ and nestedParentToChildren_.
///  4. refuseThrowsLeavingBoxes: nothing; it refuses a fault or pinch-out inside
///     a box that joins a box cell to a cell outside the box.
///  5. findDroppedFaces: faceDropped_.
///  6. findFaultedSides: faultedSides_ and acrossRange_.
///  7. buildLeafFaces: leafFaces_, mosaicOutside_, refinedCartToCell_,
///     numSourceFaces_, and each box's faceToLeaf.
///  8. buildSyntheticFaces: syntheticFaces_ and numLeafFaces_.
///  9. buildLeafTopology: the leaf's face_to_cell, face_to_point, cell_to_face,
///     face tags, normals and face geometries.
/// 10. computeCellData: the leaf's cell_to_point and cell geometries,
///     leafGlobalCell_, leafChildToParent_, leafCartesianIndexInParentCell_ and
///     leafToLevel_.
/// 11. buildParallelStructures: the leaf's parallel cell index set, partition
///     types, communication interfaces and global ids (distributed grids only).
/// 12. recordLeafMetadata: the leaf's metadata, the leaf index of the cells of
///     level zero and of the level grids, the parent-to-children relations of
///     level zero and of the parent LGRs, and the levels' refinement maximum.
///
/// A later stage modifies a member of an earlier one only here: buildCornerPool
/// and buildLeafFaces write their fields of boxes_; buildSyntheticFaces adds the
/// split corners to the leaf and updates numLeafCorners_; every stage after
/// buildCornerPool fills the leaf; and recordLeafMetadata moves the relations it
/// records out of their members. outsideNeighborOf() caches its answers.
///
/// The members of one assembly are filled once: build() runs once (asserted).
class LeafGridAssembler
{
public:
    /// Where a leaf cell or face comes from. Public, as a helper of the
    /// implementation reads the leaf faces.
    struct SourceRef
    {
        int grid;  // 0 = level zero, b+1 = level grid of box b
        int index; // cell or face index in that grid
    };

    /// Stores the input, and computes the box tree and the basic level-zero
    /// lookups.
    ///
    /// @param storage   The hierarchy [level0, level1..levelB]; level grids as
    ///                  produced by assembleBlockLevelGrid for requests[b].
    /// @param requests  The block requests, one per level grid; a parent LGR comes
    ///                  before the boxes nested in it.
    /// @param parentDims, coord, zcorn, actnum  The parent corner-point
    ///                  description (actnum may be null).
    /// @param comm      Communicator for the leaf grid object.
    LeafGridAssembler(const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& storage,
                      const std::vector<BlockRefinement>& requests,
                      const std::array<int,3>& parentDims,
                      const double* coord,
                      const double* zcorn,
                      const int* actnum,
                      Dune::MPIHelper::MPICommunicator comm);

    /// Runs the stages; call once.
    /// @return the assembled leaf; the caller appends it to storage.
    std::shared_ptr<Dune::cpgrid::CpGridData> build();

private:
    // TODO: store the box's level index (b + 1) here, and give SourceRef a boxIndex(), to
    // drop the b + 1 and grid - 1 conversions in the stages.
    struct BoxData
    {
        Dune::cpgrid::CpGridData* level;
        int cornerMask{0};      // cornerSlotMask of the level grid
        std::array<int,3> factors;
        std::array<AxisSubdivision,3> subs;
        std::array<AxisPositions,3> pos;
        bool graded{false};
        std::array<int,3> refinedDims;
        // children of each parent cell, ordered by cartesianIndexInParentCell
        std::map<int, std::vector<int>> childrenOfParent;
        // level corner -> corner of the parent grid (or -1)
        std::vector<int> cornerEquiv;
        // level corner -> leaf pool index
        std::vector<int> cornerToLeaf;      // written by buildCornerPool
        // level cell/face -> leaf index
        std::vector<int> cellToLeaf;
        std::vector<int> faceToLeaf;        // written by buildLeafFaces
    };

    struct SyntheticFace
    {
        int boxCell;
        int neighborCell;                        // coarse, or a child of another box
        std::vector<int> points;                 // leaf corner indices
        Dune::FieldVector<double,3> normal;      // unit, along +axis (CpGrid convention)
        Dune::FieldVector<double,3> center;
        double area;
        enum face_tag tag;
        int side;                                // box boundary side: -1 low, +1 high
    };

    // The constructor's lookups; each reads the members declared before the one it
    // initialises.
    std::vector<int> computeParentBoxOf() const;
    std::vector<int> computeCompressed0() const;
    std::vector<int> computeBoxOfCell() const;

    // The stages, in the order build() runs them.
    void collectBoxData();
    void buildCornerPool();
    void buildLeafCells();
    void refuseThrowsLeavingBoxes() const;
    void findDroppedFaces();
    void findFaultedSides();
    void buildLeafFaces();
    void buildSyntheticFaces();
    void buildLeafTopology();
    void computeCellData();
    void buildParallelStructures();
    void recordLeafMetadata();

    // Used by several stages.
    /// Whether box b is nested, refining another box's level grid.
    bool nested(int b) const;
    /// The box refining a cell of the grid box b refines, or -1 (after collectBoxData).
    int boxOfParentGridCell(int b, int cell) const;
    /// The leaf index of an unrefined cell of the grid box b refines (after buildLeafCells).
    int leafOfParentGridCell(int b, int cell) const;
    /// The IJK of a cell of the grid box b refines (after collectBoxData).
    std::array<int,3> parentGridIJK(int b, int cell) const;
    /// The cell across a nested box's side (after collectBoxData).
    int nestedNeighborOf(int b, int parent, int axis, int side) const;
    /// The level-zero cell across a host cell's face on a box side (after collectBoxData).
    int outsideNeighborOf(int parent, int axis, int side) const;

    // Const access for the stages after the one that writes each member.
    const std::vector<BoxData>& boxes() const { return boxes_; }
    const std::vector<std::vector<int>>& nestedBoxOf() const { return nestedBoxOf_; }
    const std::vector<SourceRef>& leafCells() const { return leafCells_; }
    const std::vector<int>& leafIdxOfCell0() const { return leafIdxOfCell0_; }
    const std::vector<std::vector<int>>& leafIdxOfLevelCell() const { return leafIdxOfLevelCell_; }
    const std::vector<std::tuple<int, std::vector<int>>>& level0ToChildren() const { return level0ToChildren_; }
    int numLeafCells() const { return numLeafCells_; }
    const std::vector<bool>& faceDropped() const { return faceDropped_; }
    const std::set<std::tuple<int,int,int>>& faultedSides() const { return faultedSides_; }
    const std::map<std::tuple<int,int,int>, CellRange>& acrossRange() const { return acrossRange_; }
    const std::vector<SourceRef>& leafFaces() const { return leafFaces_; }
    const std::vector<int>& mosaicOutside() const { return mosaicOutside_; }
    const std::vector<std::map<int,int>>& refinedCartToCell() const { return refinedCartToCell_; }
    int numSourceFaces() const { return numSourceFaces_; }
    const std::vector<SyntheticFace>& syntheticFaces() const { return syntheticFaces_; }
    int numLeafFaces() const { return numLeafFaces_; }
    int numLeafCorners() const { return numLeafCorners_; }

    // The input.
    const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& storage_;
    const std::vector<BlockRefinement>& requests_;
    const std::array<int,3> parentDims_;
    const double* const coord_;
    const double* const zcorn_;
    const int* const actnum_;
    const Dune::MPIHelper::MPICommunicator comm_;

    // Written by the constructor: the box tree and the basic level-zero lookups.
    Dune::cpgrid::CpGridData& level0_;
    const int numBoxes_;
    const std::array<int,3> dims0_;
    const int numCells0_;
    const int numCorners0_;
    // The parent box of each box, -1 for a box that refines level zero.
    const std::vector<int> parentBoxOf_;
    // The level-zero cell at each Cartesian index of level zero, or -1.
    const std::vector<int> compressed0_;
    // The box of each level-zero cell, or -1.
    const std::vector<int> boxOfCell_;
    const std::vector<std::array<int,8>>& cellToPoint0_;
    const Dune::cpgrid::OrientedEntityTable<1,0>& faceToCell0_;
    const Opm::SparseTable<int>& faceToPoint0_;
    const Dune::cpgrid::EntityVariable<enum face_tag, 1>& faceTag0_;
    const Dune::cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>, 1>& faceNormals0_;
    const int numFaces0_;

    // Written by collectBoxData.
    std::vector<BoxData> boxes_;
    int mask0_{0};      // cornerSlotMask of level zero
    // [parent box][cell of its level grid] -> nested (child) box index, or -1
    std::vector<std::vector<int>> nestedBoxOf_;

    // Written by buildCornerPool. The later stages fill the leaf.
    std::shared_ptr<Dune::cpgrid::CpGridData> leaf_;
    // leafCornerOrigin_[corner] = {level, index}: the grid where the leaf corner first
    // appears, and its index there; {0, corner} for the corners of level zero.
    std::vector<std::array<int,2>> leafCornerOrigin_;
    int numLeafCorners_{0};     // updated by buildSyntheticFaces

    // Written by buildLeafCells.
    std::vector<SourceRef> leafCells_;
    std::vector<int> leafIdxOfCell0_;
    std::vector<std::vector<int>> leafIdxOfLevelCell_;
    // level0ToChildren_[c] = {level, children}: the immediate children of level-zero cell c in
    // its box's level grid ({level, {}} if the box's MINPV removed them all), or {-1, {}}.
    // TODO: nested vectors are slow to fill and read; a SparseTable<int> with the level stored
    // alongside would do better, here and for nestedParentToChildren_.
    std::vector<std::tuple<int, std::vector<int>>> level0ToChildren_;
    int numLeafCells_{0};
    std::map<int, std::vector<std::tuple<int, std::vector<int>>>> nestedParentToChildren_;

    // Written by findDroppedFaces: by level-zero face, whether the leaf drops it.
    std::vector<bool> faceDropped_;

    // Written by outsideNeighborOf(): its answers, by (level-zero cell, axis 0..2, side -1
    // or +1); top-level boxes only.
    mutable std::map<std::tuple<int,int,int>, int> verifiedNeighbor_;

    // Written by findFaultedSides: the faulted sides, as (top-level box index, axis 0..2,
    // side -1 or +1).
    std::set<std::tuple<int,int,int>> faultedSides_;
    // By the same (box, axis, side): the parent IJK range of the cells really across that
    // box side, from the level-zero faces: a throw taller than the box puts them outside
    // the box's own layers.
    std::map<std::tuple<int,int,int>, CellRange> acrossRange_;

    // Written by buildLeafFaces.
    std::vector<SourceRef> leafFaces_;
    // The leaf cell across a face that its source grid sees as a box boundary face (an
    // unrefined cell of the grid the box refines, or a child of another box), aligned
    // with leafFaces_; -1 otherwise.
    std::vector<int> mosaicOutside_;
    std::vector<std::map<int,int>> refinedCartToCell_;
    int numSourceFaces_{0};

    // Written by buildSyntheticFaces.
    std::vector<SyntheticFace> syntheticFaces_;
    int numLeafFaces_{0};

    // Written by computeCellData; recordLeafMetadata moves them into the leaf.
    std::vector<int> leafGlobalCell_;
    std::vector<std::array<int,2>> leafChildToParent_;
    std::vector<int> leafCartesianIndexInParentCell_;
    std::vector<std::array<int,2>> leafToLevel_;

    bool built_{false};
};

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_LEAFGRIDASSEMBLER_HEADER_INCLUDED
