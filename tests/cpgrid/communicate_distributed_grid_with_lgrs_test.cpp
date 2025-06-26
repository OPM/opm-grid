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
#include "config.h"

#define BOOST_TEST_MODULE CommunicateDistributedGridWithLgrsTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

// Only for testing - copied from opm-simulators/opm/models/parallel/gridcommhandles.hh

namespace Opm {

/*!
 * \brief Data handle for parallel communication which sums up all
 *        values are attached to DOFs
 */
template <class FieldType, class Container, class EntityMapper, int commCodim>
class GridCommHandleSum
    : public Dune::CommDataHandleIF<GridCommHandleSum<FieldType, Container,
                                                      EntityMapper, commCodim>,
                                    FieldType>
{
public:
    GridCommHandleSum(Container& container, const EntityMapper& mapper)
        : mapper_(mapper), container_(container)
    {}

    bool contains(int, int codim) const
    {
        // return true if the codim is the same as the codim which we
        // are asked to communicate with.
        return codim == commCodim;
    }

    bool fixedSize(int, int) const
    {
        // for each DOF we communicate a single value which has a
        // fixed size
        return true;
    }

    template <class EntityType>
    std::size_t size(const EntityType&) const
    {
        // communicate a field type per entity
        return 1;
    }

    template <class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp& buff, const EntityType& e) const
    {
        const unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));
        buff.write(container_[dofIdx]);
    }

    template <class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp& buff, const EntityType& e, std::size_t)
    {
        const unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));

        FieldType tmp;
        buff.read(tmp);
        container_[dofIdx] += tmp;
    }

private:
    const EntityMapper& mapper_;
    Container& container_;
};
}

struct Fixture {
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

void createTestGrid(Dune::CpGrid& grid)
{
    Opm::Parser parser;
    const std::string deck_string = R"(
RUNSPEC
DIMENS
  4 4 1 /
GRID
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  2  3  2  3  1  1  6  6  3/
ENDFIN
DX
  16*1000 /
DY
	16*1000 /
DZ
	16*20 /
TOPS
	16*8325 /
ACTNUM
        16*1 /
PORO
  16*0.15 /
PERMX
  16*1 /
COPY
  PERMX PERMZ /
  PERMX PERMY /
/
EDIT
OIL
GAS
TITLE
The title
START
16 JUN 1988 /
PROPS
REGIONS
SOLUTION
SCHEDULE
)";

    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();
    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);
}

void checkChildGlobalIdsTest(const Dune::CpGrid& grid)
{
    // first_child_ids [ parent cell global id ] =  first child cell global id from undistributed view.
    const std::unordered_map<int,int> first_child_ids = { {5, 66}, {6, 93}, {9, 120}, {10, 147} };
    // parent cell with global id equal to 5 has children cell ids 66, ..., 92 (= 66 + 26).
    // parent cell with global id equal to 6 has children cell ids 93, ..., 119 (= 119 + 26).
    // parent cell with global id equal to 9 has children cell ids 120, ..., 146 (= 120 + 26).
    // parent cell with global id equal to 10 has children cell ids 147, ..., 173 (= 147 + 26).

    for (const auto& element : elements(grid.leafGridView())) {
        if (!element.hasFather()) {
            continue;
        }
        const auto& parent_globalId = grid.globalIdSet().id(element.father());
        const auto& idx_in_parent = element.getIdxInParentCell();

        const auto& expected_elem_globalId = first_child_ids.at(parent_globalId) + idx_in_parent;
        const auto& actual_elem_globalId = grid.globalIdSet().id(element);
        BOOST_CHECK_EQUAL( expected_elem_globalId, actual_elem_globalId);
    }
}

void checkCommunicationWorks(const Dune::CpGrid& grid)
{
    // This pretends to do something on the lines
    // of fvbasadiscretization finishInit()

    double gridTotalVolume = 0.0;
    std::vector<double> dofTotalVolume{};

    const auto& gridView = grid.leafGridView();
    // initialize the volume of the finite volumes to zero
    const std::size_t numDof = gridView.size(0);
    dofTotalVolume.resize(numDof, 0.0);

    for (const auto& element : Dune::elements(gridView)) {
        // ignore everything which is not in the interior if the
        // current process' piece of the grid
        if (element.partitionType() != Dune::InteriorEntity) {
            continue;
        }

        // random computation - only for testing
        const double dofVolume = element.geometry().volume();
        dofTotalVolume[element.index()] += dofVolume;
        gridTotalVolume += dofVolume;
    }

    const Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> leafMapper(gridView, Dune::mcmgElementLayout());

    // add the volumes of the DOFs on the process boundaries
    Opm::GridCommHandleSum<double,
                           std::vector<double>,
                           Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView>,
                           0> sumHandle(dofTotalVolume, leafMapper);

    gridView.communicate(sumHandle,
                         Dune::InteriorBorder_All_Interface,
                         Dune::ForwardCommunication);

    // After communication, all DOF volumes should be strictly positive.
    for (const auto& volume : dofTotalVolume) {
        BOOST_CHECK_GT(volume, 0.0); // GT: greater than
    }
}

BOOST_AUTO_TEST_CASE(callSyncCellIdsAndCommunicate)
{
    // Create the grid and add the LGRs in the global-view (non-distributed view)
    Dune::CpGrid grid;
    createTestGrid(grid);
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3, 3, 3}},
                               /* startIJK_vec = */ {{1, 1, 0}},
                               /* endIJK_vec = */  {{3, 3, 1}},
                               /* lgr_name_vec = */  {"LGR1"});

    checkChildGlobalIdsTest(grid);

    if (grid.comm().size()>1) {
        // Distribute level zero grid
        grid.loadBalance(/*overlapLayers*/ 1,
                         /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                         /*imbalanceTol*/ 1.1,
                         /*level*/ 0);

        // Add LGRs on the distributed view of level zero
        grid.addLgrsUpdateLeafView(/*cells_per_dim_vec*/ {{3, 3, 3}},
                                   /*startIJK_vec*/ {{1, 1, 0}},
                                   /*endIJK_vec*/ {{3, 3, 1}},
                                   /*lgr_name_vec*/ {"LGR1"});

        // CpGrid::syncDistributedGlobalCellIds() rewrite the cell ids of refined
        // cells to coincide with the ones from the undistributed view.
        grid.syncDistributedGlobalCellIds();
        checkChildGlobalIdsTest(grid);

        checkCommunicationWorks(grid);
    }
}

BOOST_AUTO_TEST_CASE(doNOTcallSyncCellIdsAndCommunicate)
{
    // Create the grid and add the LGRs in the global-view (non-distributed view)
    Dune::CpGrid grid;
    createTestGrid(grid);
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3, 3, 3}},
                               /* startIJK_vec = */ {{1, 1, 0}},
                               /* endIJK_vec = */  {{3, 3, 1}},
                               /* lgr_name_vec = */  {"LGR1"});

    checkChildGlobalIdsTest(grid);

    if (grid.comm().size()>1) {

        // Load balance the grid and add the LGRs in the distributed view.
        grid.loadBalance(/*overlapLayers*/ 1,
                         /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                         /*imbalanceTol*/ 1.1,
                         /*level*/ 0);
        grid.addLgrsUpdateLeafView(/*cells_per_dim_vec*/ {{3, 3, 3}},
                                   /*startIJK_vec*/ {{1, 1, 0}},
                                   /*endIJK_vec*/ {{3, 3, 1}},
                                   /*lgr_name_vec*/ {"LGR1"});

        // Commented out to test if communication occurs regardless of cell id synchronization.
        // grid.syncDistributedGlobalCellIds();
        // This fails because cell ids are not synchronized across processes.
        // checkChildGlobalIdsTest(grid);

        // Communication is possible even without synchronized cell ids across processes.
        checkCommunicationWorks(grid);
    }
}
