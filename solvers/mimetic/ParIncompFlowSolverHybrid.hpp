//===========================================================================
//
// File: ParIncompFlowSolverHybrid.hpp
//
// Created: Tue Jun 30 10:25:40 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009 SINTEF ICT, Applied Mathematics.
  Copyright 2009 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_PARINCOMPFLOWSOLVERHYBRID_HEADER
#define OPENRS_PARINCOMPFLOWSOLVERHYBRID_HEADER

#include "config.h"
#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <tr1/unordered_map>
#include <boost/bind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/SparseTable.hpp>
#include <dune/common/StopWatch.hpp>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/io.hh>
#include <dune/istl/indexset.hh>
#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <dune/solvers/common/BoundaryConditions.hpp>
#include <dune/solvers/common/Matrix.hpp>


namespace Dune {
    namespace {
        /// @brief
        ///    Verify that each pair of connected cells in a grid are
        ///    connected by a single face only.
        ///
        /// @tparam GI
        ///    Type presenting an interface to a grid (typically a
        ///    discretized geological model).  The type is assumed to
        ///    expose a forward iterator type, @code CellIter
        ///    @endcode, and a pair of delimiters @code cellbegin()
        ///    @endcode and @code cellend() @endcode to traverse the
        ///    cells of a grid.  The cell iterator, in turn, is
        ///    expected to expose a forward iterator type @code
        ///    FaceIter @endcode, and a pair of delimiters @code
        ///    facebegin() @endcode and @code faceend() @endcode to
        ///    traverse the faces of a cell.
        ///
        /// @param [in] g
        ///    The grid.
        ///
        /// @return
        ///    @code true @endcode if every pair of connected grid
        ///    cells are connected by a single face only and @code
        ///    false @endcdode otherwise.
        template<class GI>
        bool topologyIsSane(const GI& g)
        {
            typedef typename GI::CellIterator CI;
            typedef typename CI::FaceIterator FI;

            bool sane = g.numberOfCells() >= 0;

            for (CI c = g.cellbegin(); sane && c != g.cellend(); ++c) {
                std::vector<int> n;

                for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                    if (!f->boundary()) {
                        n.push_back(f->neighbourCellIndex());
                    }
                }
                std::sort(n.begin(), n.end());

                sane = std::unique(n.begin(), n.end()) == n.end();
            }

            return sane;
        }


        /// @brief
        ///    Type for constructing a binary function object which,
        ///    together with the standard transforming algorithm @code
        ///    std::transform() @endcode, implements the extended BLAS
        ///    level 1 operation @f[ y \leftarrow ax + by. @f]
        ///
        /// @details
        ///    Suppose @f$x@f$ and @f$y@f$ are collections of the same
        ///    @code size() @endcode, for instance two equally sized
        ///    @code std::vector<T>@endcode's for some element type
        ///    @code T @endcode.  Then the statements
        ///    @code
        ///       T a = <some value>;
        ///       T b = <some other value>;
        ///       std::transform(x.begin(), x.end(), y.begin(),
        ///                      y.begin(), axpby<T>(a, b));
        ///    @endcode
        ///    implement the extended BLAS level 1 operation outlined
        ///    above for the elements of @f$x@f$ and @f$y@f$.
        ///
        /// @tparam T
        ///    Element type of the collections to which the functor
        ///    will be applied.  Expected to be an arithmetic type,
        ///    and usually @code T @endcode is an alias for @code
        ///    double @endcode.
        template<typename T>
        class axpby : public std::binary_function<T,T,T> {
        public:
            /// @brief Constructor.
            /// @param [in] a
            ///    Constant multiplying the elements of @f$x@f$.
            /// @param [in] b
            ///    Constant multiplying the elements of @f$y@f$.
            axpby(const T& a, const T& b) : a_(a), b_(b) {}

            /// @brief
            ///    Function call operator required for functor
            ///    behaviour.
            ///
            /// @param [in] x
            ///    Single element of collection @f$x@f$.
            ///
            /// @param [in] y
            ///    Single element of collection @f$y@f$.
            ///
            /// @return
            ///    Corresponding transformed value @f$ax + by@f$.
            T operator()(const T& x, const T& y)
            {
                return a_*x + b_*y;
            }
        private:
            T a_, b_;
        };
    }


    /// @brief
    ///    Solve mixed formulation of incompressible flow modelled by
    ///    Darcy's law
    ///    @f[@f{aligned}{
    ///       v &= -K\lamda(\nabla p + \rho\vec{g}), \\ \nabla\cdot v &= q.
    ///    @f}@f]
    ///    The solver is based on a hybrid discretization scheme which
    ///    yields a system of linear equations of the form
    ///    @f[@f{bmatrix}{
    ///        B & C & D \\ C^{T} & 0 & 0 \\ D^{T} & 0 & 0
    ///      }\,
    ///      @f{bmatrix}{v \\ -p \\ \pi@f} =
    ///      @f{bmatrix}{f \\ g \\ h@f}
    ///    @f]
    ///    where @f$v@f$ represents the interface fluxes for each
    ///    interface for each cell, @f$p@f$ are the cell pressures and
    ///    @f$\pi@f$ are the interface pressures.
    ///
    ///    Through a Schur complement analysis, the above system of
    ///    linear equations is reduced to an equivalent system of the
    ///    form
    ///    @f[@f{bmatrix}{
    ///        B & C & D \\ 0 & -L & -F \\ 0 & 0 & S
    ///      }\,
    ///      @f{bamtrix}{v \\ -p \\ \pi@f} =
    ///      @f{bmatrix}{f \\ \hat{g} \\ r@f}
    ///    @f]
    ///    where @f$L = C^{T}B^{-1}C@f$, @f$F=C^{T}B^{-1}D@f$, and
    ///    @f$S = D^{T}B^{-1}D - F^{T}L^{-1}F@f$.  Similarly,
    ///    @f$\hat{g} = g - C^{T}B^{-1}f@f$ and @f$r = D^{T}B^{-1}f +
    ///    F^{T}L^{-1}\hat{g} - h@f$.
    ///
    ///    The system @f$S\pi = r@f$ is the system of linear equations
    ///    which is actually solved using separte linear system solver
    ///    software.  Then, a back substitution process yields the
    ///    cell pressures and interface fluxes by solving the simpler
    ///    systems
    ///    @f[@f{aligned}{
    ///       Lp &= \hat{g} + F\pi \\ Bv &= f + Cp - D\pi
    ///    @f}@f]
    ///    Specifically, the matrix @f$L@f$ is diagonal (a single
    ///    scalar per grid cell) and the matrix @f$B^{-1}@f$ is
    ///    available as a bi-product of the Schur reduction process.
    ///    Finally, in the case of mimetic discretizations using the
    ///    @code MimeticIPEvaluator @endcode class, the matrix
    ///    @f$B^{-1}@f$ is directly available through an explicit
    ///    formula.
    ///
    ///    The matrix @f$B@f$ is a product of two parts--one which
    ///    depends only on the grid (i.e., the topology and geomtry)
    ///    and the geological properties (i.e., the permeability
    ///    field) and one part which depends only on the temporally
    ///    varying fluid data (i.e., viscosities and saturation
    ///    dependent relative permeabilities).  Furthermore, the
    ///    matrices @f$C@f$ and @f$D@f$ depend @em only on the grid
    ///    topolgy.  Consequently, the cost of assembling the system
    ///    @f$S@f$ may be reduced, at least for simulations involving
    ///    multiple pressure solves, by performing all of the static
    ///    (i.e., the parts dependent only on the grid and geological
    ///    properties) assembly once, at system initalization.  Then
    ///    the system assembly process consists of calculating the
    ///    saturation dependent contributions and adding these into
    ///    the system of linear equations.
    ///
    ///    The code is structured in a manner similar to traditional
    ///    FEM software implementations.  In particular, we assemble
    ///    the coefficient matrix @f$S@f$ and system right hand side
    ///    @f$r@f$ on a cell-by-cell basis.  This yields a number of
    ///    important simplifications.  In particular, the matrix
    ///    @f$D@f$ reduces to the identity on a single cell.  Its
    ///    effects are produced only through a local-to-global map of
    ///    degrees of freedom.  Similarly, on a single grid cell, the
    ///    matrix @f$C=[1,1,\dots,1]^T@f$ with the number of elements
    ///    equal to the number of interfaces (i.e., neighbours) of the
    ///    cell.  Finally, the matrix @f$B^{-1}@f$ is computed as
    ///    @f[B^{-1} = \hat{B}^{-1}/\lambda_T@f] where
    ///    @f$\hat{B}^{-1}@f$ denotes the static part of @f$B^{-1}@f$
    ///    and @f$\lambda_T@f$ denotes the @em total fluid @em
    ///    mobility in the grid cell.  For ease of implementation of
    ///    the back substitution process, we compute and store the
    ///    values of @f$L@f$, @f$F@f$ and @f$\hat{g}@f$ during each
    ///    system assembly process.
    ///
    ///    A final quirk of the implementation is that we always solve
    ///    for every interface pressure, even if a given pressure
    ///    value is known through, e.g., a prescribed pressure
    ///    boundary condition.  This feature enables changing the type
    ///    and value of a set of boundary conditions between pressure
    ///    solves without having to reconstruct the sparsity structure
    ///    of the coefficient matrix @f$S@f$--at least while no
    ///    connections are introduced or removed.  Periodic boundary
    ///    condtions introduce new connections (i.e., new non-zero
    ///    coefficient matrix entries) so in order to gain maximum
    ///    efficiency, you should initialize the flow solver system
    ///    only after you know the existence of all periodic boundary
    ///    conditions which will be employed in a given simulation.
    ///
    ///    The use of non-periodic boundary conditions after periodic
    ///    connections have been introduced is fully supported.  This
    ///    mode introduces a performance caveat.  Specifically, a few
    ///    entries which would otherwise not have been entered into
    ///    the matrix due to being automatically zero will become
    ///    explicit zeros instead.  Consequently, each iteration of an
    ///    iterative solver for the system @f$S\pi = r@f$ becomes
    ///    slightly more expensive as explicit zero-multiplications
    ///    occur.  Moreover, introducing explicit non-zero
    ///    representation of matrix entries which are always zero
    ///    increases the demand for memory resources.
    ///
    ///    The flow solver is initialized in a three-step process.
    ///    The first step, represented by method @code clear()
    ///    @endcode, releases any previously held data in the internal
    ///    data structures and prepares the solver system for defining
    ///    a new problem.  The second step, represented by the @code
    ///    initSystemStructure() @endcode enumerates the primary
    ///    degrees of freedom (i.e., the interface pressure values)
    ///    and determines, and allocates, the coefficient matrix
    ///    non-zero structure.  The third and final step, represented
    ///    by method @code computeInnerProducts() @endcode, computes
    ///    the static (i.e., geology and geomtry-dependent) inner
    ///    product matrices @f$\hat{B}^{-1}@f$.  Method @code
    ///    computeInnerProducts() @endcode is offered separately in
    ///    order to support solve several different property models on
    ///    the same underlying geometry (grid).  Finally, method @code
    ///    init() @endcode is a convenience method which calls the
    ///    other initialization methods in sequence.
    ///
    ///    Following solver intialization, a sequence of flow problems
    ///    differing only by boundary condition type/value and/or
    ///    differing fluid saturation values may be resolved by
    ///    separate calls to the @code solve() @endcode method.  At
    ///    any time following a call to @code solve() @endcode may the
    ///    solution, represented by the type @code SolutionType
    ///    @endcode, to the most recently resolved flow problem be
    ///    retrieved through the @code getSolution() @endcode method.
    ///    We note that @code SolutionType @endcode is @em always a
    ///    reference-to-const.
    ///
    /// @tparam GridInterface
    ///    Type presenting an interface to a grid (typically a
    ///    discretized geological model).  The type is assumed to
    ///    expose a forward iterator type, @code CellIter @endcode,
    ///    and a pair of delimiters @code cellbegin() @endcode and
    ///    @code cellend() @endcode to traverse the cells of a grid.
    ///    The cell iterator, in turn, is expected to expose a forward
    ///    iterator type @code FaceIter @endcode, and a pair of
    ///    delimiters @code facebegin() @endcode and @code faceend()
    ///    @endcode to traverse the faces of a cell.
    ///
    /// @tparam ReservoirInterface
    ///    Type presenting an interface to reservoir properties such
    ///    as permeability, porosity, and various fluid
    ///    characteristics (density, mobility &c).  The type is
    ///    assumed to expose a method, @code permeability() @endcode
    ///    through which the assigned permeability of a single grid
    ///    cell, represented by a @code GridInterface::CellIter
    ///    @endcode, may be recovered.  The type is further expected
    ///    to provide methods @code phaseMobility() @endcode and @code
    ///    phaseDensity() @endcode for phase mobility and density in a
    ///    single cell, respectively.
    ///
    /// @tparam BCInterface
    ///    Type presenting an interface to boundary conditions.  The
    ///    type is expected to provide a method, @code flowCond()
    ///    @endcode, from which a @code BCInterface::FlowBC @endcode
    ///    boundary condtion type of any given boundary face may be
    ///    recovered.  The flow boundary condition type is expected to
    ///    provide methods for quering the type of boundary condition
    ///    (i.e., prescribed pressure values, prescribed flux values,
    ///    or prescribed pressure drops in the case of periodic
    ///    boundary conditions) as well as the numerical values of
    ///    these conditions.
    ///
    /// @tparam InnerProduct
    ///    Type presenting a specific inner product defining a
    ///    discretization of the Darcy equation.
    template<class                          GridInterface,
             class                          ReservoirInterface,
             class                          BCInterface,
             template<class,int,bool> class InnerProduct>
    class ParIncompFlowSolverHybrid {
        /// @brief
        ///    The element type of the matrix representation of the
        ///    mimetic inner product.  Assumed to be a floating point
        ///    type, and usually, @code Scalar @endcode is an alias
        ///    for @code double @endcode.
        typedef typename GridInterface::Scalar Scalar;

        /// @brief
        ///   Tag for identifying which type of condition/equation to
        ///   assemble into the global system of linear equation for
        ///   any given face/connection.
        enum FaceType { Internal, Dirichlet, Neumann, Periodic };

        /// @brief
        ///    Type representing the solution to a given flow problem.
        class FlowSolution {
        public:
            /// @brief
            ///    The element type of the matrix representation of
            ///    the mimetic inner product.  Assumed to be a
            ///    floating point type, and usually, @code Scalar
            ///    @endcode is an alias for @code double @endcode.
            typedef typename GridInterface::Scalar       Scalar;

            /// @brief
            ///    Convenience alias for the grid interface's cell
            ///    iterator.
            typedef typename GridInterface::CellIterator CI;

            /// @brief
            ///    Convenience alias for the cell's face iterator.
            typedef typename CI           ::FaceIterator FI;

            friend class ParIncompFlowSolverHybrid;

            /// @brief
            ///    Retrieve the current cell pressure in a given cell.
            ///
            /// @param [in] c
            ///    Cell for which to retrieve the current cell
            ///    pressure.
            ///
            /// @return
            ///    Current cell pressure in cell @code *c @endcode.
            Scalar pressure(const CI& c) const
            {
                return pressure_[cellno_[c->index()]];
            }

            /// @brief
            ///    Retrieve current flux across given face in
            ///    direction of outward normal vector.
            ///
            /// @param [in] f
            ///    Face across which to retrieve the current outward
            ///    flux.
            ///
            /// @return
            ///    Current outward flux across face @code *f @endcode.
            Scalar outflux (const FI& f) const
            {
                return outflux_[cellno_[f->cellIndex()]][f->localIndex()];
            }
        private:
            std::vector< int  > cellno_;
            SparseTable< int  > cellFaces_;
            std::vector<Scalar> pressure_;
            SparseTable<Scalar> outflux_;
        };

    public:

        /// @brief Default constructor.
        ParIncompFlowSolverHybrid()
	    : comm_(MPI_COMM_WORLD)
        {
            clear();
        }


        /// @brief
        ///    All-in-one initialization routine.  Enumerates all grid
        ///    connections, allocates sufficient space, defines the
        ///    structure of the global system of linear equations for
        ///    the contact pressures, and computes the permeability
        ///    dependent inner products for all of the grid's cells.
        ///
        /// @param [in] g
        ///    The grid.
        ///
        /// @param [in] r
        ///    The reservoir properties of each grid cell.
        ///
        /// @param [in] bc
        ///    The boundary conditions describing how the current flow
        ///    problem interacts with the outside world.  This is used
        ///    only for the purpose of introducing additional
        ///    couplings in the case of periodic boundary conditions.
        ///    The specific values of the boundary conditions are not
        ///    inspected in @code init() @endcode.
        void init(const GridInterface&      g,
                  const ReservoirInterface& r,
                  const BCInterface&        bc,
                  const std::vector<int>&   partition,
		  const int my_partition)
        {
            clock_.start();
            clear();

            if (g.numberOfCells() > 0) {
                initSystemStructure(g, bc, partition, my_partition);
                computeInnerProducts(r);
                printElapsedTime("computeInnerProducts()");
            }
        }


        /// @brief
        ///    Clear all topologic, geometric and rock-dependent
        ///    information currently held in internal data structures.
        ///    Release all memory resources in preparation of
        ///    constructing a solver for a new problem.  Method @code
        ///    clear() @endcode must be called prior to any other
        ///    method of the class.
        void clear()
        {
            pgrid_                  =  0;
            ppartition_             =  0;
	    my_partition_           = -1;
            max_ncf_                = -1;
            total_num_faces_        =  0;
            total_num_cells_        =  0;
            matrix_structure_valid_ = false;
            do_regularization_      = true; // Assume pure Neumann by default.

            bdry_id_map_.clear();
            ppartner_dof_.clear();

            std::vector<Scalar>().swap(L_);
            std::vector<Scalar>().swap(g_);
            Binv_.clear();   F_.clear();   f_.clear();

            std::vector<int>().swap(flowSolution_.cellno_);
            flowSolution_.cellFaces_.clear();

            cleared_state_ = true;
        }


        /// @brief
        ///    Compute structure of coefficient matrix in final system
        ///    of linear equations for this flow problem.  Allocates
        ///    all storage needed by the flow solver itself.
        ///
        /// @param [in] g
        ///    The grid.
        ///
        /// @param [in] bc
        ///    The boundary conditions describing how the current flow
        ///    problem interacts with the outside world.  This is used
        ///    only for the purpose of introducing additional
        ///    couplings in the case of periodic boundary conditions.
        ///    The specific values of the boundary conditions are not
        ///    inspected in @code init() @endcode.
        void initSystemStructure(const GridInterface& g,
                                 const BCInterface&   bc,
                                 const std::vector<int>& partition,
				 const int my_partition)
        {
            ASSERT2 (cleared_state_,
                     "You must call clear() prior to initSystemStructure()");
            ASSERT  (topologyIsSane(g));

	    my_partition_ = my_partition;
            enumerateDof(g, bc, partition);
            printElapsedTime("enumerateDof()");
            allocateConnections(bc);
            printElapsedTime("allocateConnections()");
            setConnections(bc);
            printElapsedTime("setConnections()");
        }


        /// @brief
        ///    Compute static (i.e., independent of saturation) parts
        ///    of the spatially varying inner products @f$ (v,K^{-1}w)
        ///    @f$ for each grid cell.  This is often a fairly
        ///    expensive process, so in order to increase the speed of
        ///    the system assembly we do not wish to perform this
        ///    computation for each time step unless we have to.
        ///    Moreover, if the permeability field changes but the
        ///    grid connections remain static, we do not have to build
        ///    a new matrix structure but need only compute new values
        ///    for these static inner products.
        ///
        /// @param [in] r
        ///    The reservoir properties of each grid cell.  In method
        ///    @code computeInnerProducts() @endcode, we only inspect
        ///    the permeability field of @code r @endcode.
        void computeInnerProducts(const ReservoirInterface& r)
        {
            ASSERT2 (matrix_structure_valid_,
                     "You must call connectionsCompleted() prior "
                     "to computeInnerProducts()");

            typedef typename GridInterface     ::CellIterator               CI;
            typedef typename ReservoirInterface::PermTensor                 PermTensor;
            typedef          InnerProduct<CI,GridInterface::Dimension,true> IP;

            IP ip(max_ncf_);
            int i = 0;
            const SparseTable<int>& cellFaces = flowSolution_.cellFaces_;
            for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
		if ((*ppartition_)[c->index()] != my_partition_) continue;
                const int nf = cellFaces[i].size();

                SharedFortranMatrix Binv(nf, nf, &Binv_[i][0]);

                ip.evaluate(c, r.permeability(c->index()), Binv);
		++i;
            }
        }


        /// @brief
        ///    Construct and solve system of linear equations for the
        ///    pressure values on each interface/contact between
        ///    neighbouring grid cells.  Recover cell pressure and
        ///    interface fluxes.  Following a call to @code solve()
        ///    @encode, you may recover the flow solution from the
        ///    @code getSolution() @endcode method.
        ///
        /// @param [in] r
        ///    The reservoir properties of each grid cell.  In method
        ///    @code solve() @endcode we query this object for the
        ///    phase mobilities (i.e., @code r.phaseMobility()
        ///    @endcode) and the phase densities (i.e., @code
        ///    phaseDensity() @encode) of each phase.
        ///
        /// @param [in] sat
        ///    Saturation of primary phase.  One scalar value for each
        ///    grid cell.  This parameter currently limits @code
        ///    IncompFlowSolverHybrid @endcode to two-phase flow
        ///    problems.
        ///
        /// @param [in] bc
        ///    The boundary conditions describing how the current flow
        ///    problem interacts with the outside world.  Method @code
        ///    solve() @endcode inspects the actual values of the
        ///    boundary conditions whilst forming the system of linear
        ///    equations.
        ///
        ///    Specifically, the @code bc.flowCond(bid) @endcode
        ///    method is expected to yield a valid @code FlowBC
        ///    @endcode object for which the methods @code pressure()
        ///    @endcode, @code pressureDifference() @endcode, and
        ///    @code outflux() @endcode yield valid responses
        ///    depending on the type of the object.
        ///
        /// @param [in] src
        ///    Explicit source terms.  One scalar value for each grid
        ///    cell representing the rate (in units of m^3/s) of fluid
        ///    being injected into (>0) or extracted from (<0) a given
        ///    grid cell.
        ///
        /// @param [in] grav
        ///    Gravity vector.  Its Euclidian two-norm value
        ///    represents the strength of the gravity field (in units
        ///    of m/s^2) while its direction is the direction of
        ///    gravity in the current model.  The gravity may, in
        ///    principle, be changed between calls of @code solve()
        ///    @endcode.
        ///
        /// @param [in] residual_tolerance
        ///    Control parameter for iterative linear solver software.
        ///    The iteration process is terminated when the norm of
        ///    the linear system residual is less than @code
        ///    residual_tolerance @endcode times the initial residual.
        ///
        /// @param [in] linsolver_verbosity
        ///    Level of diagnostic output from linear system solver
        ///    output.  Semantics and possible values strongly tied to
        ///    the specifics of the currently employed linear system
        ///    solver software.
        void solve(const ReservoirInterface&  r  ,
                   const std::vector<double>& sat,
                   const BCInterface&         bc ,
                   const std::vector<double>& src,
                   const typename GridInterface::CellIterator::Vector& grav,
                   double residual_tolerance = 1e-8,
                   int linsolver_verbosity = 1)
        {
            assembleDynamic(r, sat, bc, src, grav);
            printElapsedTime("assembleDynamic()");

            // printSystem("linsys_mimetic");
#if 0
            solveLinearSystem(residual_tolerance, linsolver_verbosity);
#else
            solveLinearSystemAMG(residual_tolerance, linsolver_verbosity);
#endif
            printElapsedTime("linear solver");
            computePressureAndFluxes(r, sat);
            printElapsedTime("computePressureAndFluxes()");
        }


        /// @brief
        ///    Type representing the solution to the problem defined
        ///    by the parameters to @code solve() @endcode.  Always a
        ///    reference-to-const.  The @code SolutionType @endcode
        ///    exposes methods @code pressure(c) @endcode and @code
        ///    outflux(f) @endcode from which the cell pressure in
        ///    cell @code *c @endcode and outward-pointing flux across
        ///    interface @code *f @endcode may be recovered.
        typedef const FlowSolution& SolutionType;

        /// @brief
        ///    Recover the solution to the problem defined by the
        ///    parameters to method @code solve() @endcode.  This
        ///    solution is meaningless without a previous call to
        ///    method @code solve() @endcode.
        ///
        /// @return
        ///    The current solution.
        SolutionType getSolution()
        {
            return flowSolution_;
        }


        /// @brief
        ///    Print statistics about the connections in the current
        ///    model.  This is mostly for debugging purposes and
        ///    should only rarely be used in client code.
        ///
        /// @tparam charT
        ///    Character type of output stream.
        ///
        /// @tparam traits
        ///    Character traits of @code charT @endcode.
        ///
        /// @param os
        ///    Output stream into which the statistics will be
        ///    printed.
        template<typename charT, class traits>
        void printStats(std::basic_ostream<charT,traits>& os)
        {
            os << "IncompFlowSolverHybrid<>:\n"
               << "\tMaximum number of cell faces = " << max_ncf_ << '\n'
               << "\tTotal number of faces        = " << total_num_faces_ << '\n';

            const std::vector<int>& cell = flowSolution_.cellno_;
            os << "cell index map = [";
            std::copy(cell.begin(), cell.end(),
                      std::ostream_iterator<int>(os, " "));
            os << "\b]\n";

            const SparseTable<int>& cf = flowSolution_.cellFaces_;
            os << "cell faces     =\n";
            for (int i = 0; i < cf.size(); ++i)
            {
                os << "\t[" << i << "] -> [";
                std::copy(cf[i].begin(), cf[i].end(),
                          std::ostream_iterator<int>(os, ","));
                os << "\b]\n";
            }
        }


        //template<typename charT, class traits>
        void printElapsedTime(const std::string& postfix = "")//, std::basic_ostream<charT,traits>& os = std::cout)
        {
            //os << prefix << "   Time elapsed: " << clock_.secsSinceLast() << std::endl;
            std::cout << "Time elapsed: " << clock_.secsSinceLast() << "  [in " << postfix << "]" << std::endl;
        }

        /// @brief
        ///    Output the current (static) inner products.  This is
        ///    only meaningfull following a call to method @code
        ///    computeInnerProducts() @endcode and is mostly useful
        ///    for debugging.  Only rarely should this method be used
        ///    from client code.
        ///
        /// @tparam charT
        ///    Character type of output stream.
        ///
        /// @tparam traits
        ///    Character traits of @code charT @endcode.
        ///
        /// @param os
        ///    Output stream into which the inner products will be
        ///    printed.
        template<class charT, class traits>
        void printIP(std::basic_ostream<charT,traits>& os)
        {
            const SparseTable<int>& cf = flowSolution_.cellFaces_;
            // Loop grid whilst building (and outputing) the inverse IP matrix.
            for (int c = 0; c != cf.size(); ++c) {
                const int nf = cf[c].size();
                ImmutableFortranMatrix Binv(nf, nf, &Binv_[c][0]);

                os << c << " -> Binv = [\n" << Binv << "]\n";
            }
        }


        /// @brief
        ///    Output current system of linear equations to permanent
        ///    storage in files.  One file for the coefficient matrix
        ///    and one file for the right hand side.  This is mostly
        ///    useful whilst debugging the solver and should rarely be
        ///    used from client code.
        ///
        /// @details
        ///    The system is stored in a format which is suitable for
        ///    importing into MATLAB or MATLAB-compatible software.
        ///    In particular, using the MATLAB @code LOAD @endcode and
        ///    @code SPCONVERT @endcode functions makes it easy to
        ///    reconstruct the system of linear equations from within
        ///    MATLAB.
        ///
        /// @param [in] prefix
        ///    Prefix from which file names for the coefficient matrix
        ///    and right hand side data are derived.  Specifically,
        ///    the matrix data is output to the file @code prefix +
        ///    "-mat.dat" @endcode while the right hand side data is
        ///    output to the file @code prefix + "-rhs.dat" @endcode.
        void printSystem(const std::string& prefix)
        {
            writeMatrixToMatlab(S_, prefix + "-mat.dat");

            std::string rhsfile(prefix + "-rhs.dat");
            std::ofstream rhs(rhsfile.c_str());
            std::copy(rhs_.begin(), rhs_.end(),
                      std::ostream_iterator<VectorBlockType>(rhs, "\n"));
        }

    private:
        typedef std::pair<int,int>                 DofID;
        typedef std::tr1::unordered_map<int,DofID> BdryIdMapType;
        typedef BdryIdMapType::const_iterator      BdryIdMapIterator;

        const GridInterface* pgrid_;
        const std::vector<int>* ppartition_;
        BdryIdMapType        bdry_id_map_;
        std::vector<int>     ppartner_dof_;

        time::StopWatch clock_;

        // ----------------------------------------------------------------
        bool cleared_state_;
        int  max_ncf_;
	int  total_num_cells_;
        int  total_num_faces_;
	int  my_partition_;

	typedef int LocalId;
	typedef int GlobalId;
	typedef Dune::OwnerOverlapCopyCommunication<LocalId,GlobalId> Communication;
	Communication comm_;

        // ----------------------------------------------------------------
        std::vector<Scalar> L_, g_;
        SparseTable<Scalar> Binv_, F_, f_;

        // ----------------------------------------------------------------
        // Actual, assembled system of linear equations
        typedef FieldVector<Scalar, 1   > VectorBlockType;
        typedef FieldMatrix<Scalar, 1, 1> MatrixBlockType;

        BCRSMatrix <MatrixBlockType>      S_;    // System matrix
        BlockVector<VectorBlockType>      rhs_;  // System RHS
        BlockVector<VectorBlockType>      soln_; // System solution (contact pressure)
        bool                              matrix_structure_valid_;
        bool                              do_regularization_;

        // ----------------------------------------------------------------
        // Physical quantities (derived)
        FlowSolution flowSolution_;


        // ----------------------------------------------------------------
        void enumerateDof(const GridInterface& g,
                          const BCInterface& bc,
                          const std::vector<int>& partition)
        // ----------------------------------------------------------------
        {
            enumerateGridDof(g, partition);
            enumerateBCDof(g, bc, partition);

            pgrid_ = &g;
            ppartition_ = &partition;
            cleared_state_ = false;
        }

        // ----------------------------------------------------------------
        void enumerateGridDof(const GridInterface& g,
                              const std::vector<int>& partition)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            // Initialize cell index -> iterator order mapping.
            const int nc_all = g.numberOfCells();
            std::vector<int>(nc_all, -1).swap(flowSolution_.cellno_);
            std::vector<int>& cell = flowSolution_.cellno_;

            // First pass: count things, and build dof<->faceindex maps.
	    std::vector<int> faceindex_to_dof(g.numberOfFaces(), -1);
	    std::vector<int> dof_to_faceindex;
	    std::vector<int> part_bdy;
            int cellno = 0;
            int tot_ncf = 0;  // Number of half-faces. Or sum of neighbourhood sizes.
	    int tot_ncf2 = 0; // Sum of squared neighbourhood sizes.
            for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
                const int c0 = c->index();
		if (partition[c0] != my_partition_) continue;
                ASSERT((0 <= c0) && (c0 < nc_all) && (cell[c0] == -1));
                cell[c0] = cellno;
		const typename GridInterface::Indices& cell_faces = g.faceIndices(c0);
		int ncf = cell_faces.size();
                max_ncf_  = std::max(max_ncf_, ncf);
                tot_ncf  += ncf;
                tot_ncf2 += ncf * ncf;
		for (FI f = c->facebegin(); f != c->faceend(); ++f) {
		    int fi = f->index();
		    if (faceindex_to_dof[fi] == -1) {
			faceindex_to_dof[fi] = dof_to_faceindex.size();
			dof_to_faceindex.push_back(fi);
			bool on_part_bdy = false;
			if (!f->boundary()) {
			    on_part_bdy = (partition[f->neighbourCellIndex()] != my_partition_);
			}
			part_bdy.push_back(on_part_bdy);
		    } else {
			if (dof_to_faceindex[faceindex_to_dof[fi]] != fi) {
			    THROW("Error in creating dof <-> faceindex maps.");
			}
		    }
		}
		++cellno;
            }
            ASSERT(cellno <= nc_all);
	    total_num_cells_ = cellno;  // Size of my partition.

	    // Second pass: build index set.
	    typedef Dune::OwnerOverlapCopyAttributeSet GridAttributes;
	    typedef GridAttributes::AttributeSet GridFlag;
	    typedef Dune::ParallelLocalIndex<GridFlag> LocalIndex;
	    int num_dof = dof_to_faceindex.size();
	    Communication::PIS& pis = comm_.indexSet();
	    pis.beginResize();
	    for (int i = 0; i < num_dof; ++i) {
		GridFlag flag = part_bdy[i] ? GridAttributes::overlap : GridAttributes::owner;
		// if (my_partition_ == 0) flag = GridAttributes::owner;
		bool is_public = part_bdy[i] ? true : false;
		pis.add(dof_to_faceindex[i], LocalIndex(i, flag, is_public));
	    }
	    pis.endResize();
	    comm_.remoteIndices().template rebuild<false>();


            // Third pass: build cell-to-face mapping and the structures of
	    // the Binv_, F_, f_ and outflux objects.
            flowSolution_.cellFaces_.reserve(cellno, tot_ncf);
            flowSolution_.outflux_  .reserve(cellno, tot_ncf);
            F_   .reserve(cellno, tot_ncf );
            f_   .reserve(cellno, tot_ncf );
            Binv_.reserve(cellno, tot_ncf2);
            typedef std::vector<int>::iterator VII;
            std::vector<Scalar> zeroes(max_ncf_*max_ncf_, 0.0);
            std::vector<int> cell_dofs(max_ncf_);
            for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
                const int c0 = c->index();
		if (partition[c0] != my_partition_) continue;
                ASSERT ((0 <=      c0 ) && (     c0  < nc_all) &&
                        (0 <= cell[c0]) && (cell[c0] < total_num_cells_));
                const int ncf = g.faceIndices(c0).size();
		cell_dofs.resize(ncf);
		for (int ll = 0; ll < ncf; ++ll) {
		    cell_dofs[ll] = faceindex_to_dof[g.faceIndices(c0)[ll]];
		}
		flowSolution_.cellFaces_.appendRow(cell_dofs.begin(), cell_dofs.end());
                flowSolution_.outflux_  .appendRow(zeroes.begin(), zeroes.begin() + ncf);
                F_                      .appendRow(zeroes.begin(), zeroes.begin() + ncf);
                f_                      .appendRow(zeroes.begin(), zeroes.begin() + ncf);
                Binv_                   .appendRow(zeroes.begin(), zeroes.begin() + ncf*ncf);
            }
	    total_num_faces_ = dof_to_faceindex.size();
        }


        // ----------------------------------------------------------------
        void enumerateBCDof(const GridInterface& g,
                            const BCInterface& bc,
                            const std::vector<int>& partition)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            bdry_id_map_.clear();
            for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
		if (partition[c->index()] != my_partition_) continue;
                for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                    if (f->boundary()) {
                        const int bid = f->boundaryId();
                        if (bc.flowCond(bid).isPeriodic()) {
                            DofID dof(cell[c->index()], f->localIndex());
                            bdry_id_map_.insert(std::make_pair(bid, dof));
                        }
                    }
                }
            }

            ppartner_dof_.clear();
            if (!bdry_id_map_.empty()) {
                ppartner_dof_.assign(total_num_faces_, -1);
                for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
		    if (partition[c->index()] != my_partition_) continue;
                    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                        if (f->boundary()) {
                            const int bid = f->boundaryId();
                            if (bc.flowCond(bid).isPeriodic()) {
                                const int dof1 = cf[cell[c->index()]][f->localIndex()];

                                BdryIdMapIterator j =
                                    bdry_id_map_.find(bc.getPeriodicPartner(bid));
                                ASSERT (j != bdry_id_map_.end());
                                const int dof2 = cf[j->second.first][j->second.second];

                                ppartner_dof_[dof1] = dof2;
                                ppartner_dof_[dof2] = dof1;
                            }
                        }
                    }
                }
            }
        }



        // ----------------------------------------------------------------
        void allocateConnections(const BCInterface& bc)
        // ----------------------------------------------------------------
        {
            ASSERT2 (!cleared_state_,
                     "You must call enumerateDof() prior "
                     "to allocateConnections()");
            ASSERT  (!matrix_structure_valid_);

            // Clear any residual data, prepare for assembling structure.
            S_.setSize(total_num_faces_, total_num_faces_);
            S_.setBuildMode(BCRSMatrix<MatrixBlockType>::random);

            // Compute row sizes
            for (int f = 0; f < total_num_faces_; ++f) {
                S_.setrowsize(f, 1);
            }

            allocateGridConnections();
            allocateBCConnections(bc);

            S_.endrowsizes();

            rhs_ .resize(total_num_faces_);
            soln_.resize(total_num_faces_);
        }


        // ----------------------------------------------------------------
        void allocateGridConnections()
        // ----------------------------------------------------------------
        {
            const   SparseTable<int>& cf = flowSolution_.cellFaces_;
            typedef SparseTable<int>::row_type::const_iterator fi;

            for (int c = 0; c < cf.size(); ++c) {
                const int nf = cf[c].size();
                fi fb = cf[c].begin(), fe = cf[c].end();

                for (fi f = fb; f != fe; ++f) {
                    S_.incrementrowsize(*f, nf - 1);
                }
            }
        }


        // ----------------------------------------------------------------
        void allocateBCConnections(const BCInterface& bc)
        // ----------------------------------------------------------------
        {
            // Include system connections for periodic boundary
            // conditions, if any.  We make an arbitrary choice in
            // that the face/degree-of-freedom with the minimum
            // numerical id of the two periodic partners represents
            // the coupling.  Suppose <i_p> is this minimum dof-id.
            // We then need to introduce a *symmetric* coupling to
            // <i_p> to each of the dof's of the cell *NOT* containing
            // <i_p>.  This choice is implemented in the following
            // loop by introducing couplings only when dof1 (self) is
            // less than dof2 (other).
            //
            // See also: setBCConnections() and addCellContrib().
            //
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            if (!bdry_id_map_.empty()) {
                // At least one periodic BC.  Allocate corresponding
                // connections.
                for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
		    if ((*ppartition_)[c->index()] != my_partition_) continue;
                    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                        if (f->boundary()) {
                            const int bid = f->boundaryId();
                            if (bc.flowCond(bid).isPeriodic()) {
                                // dof-id of self
                                const int dof1 = cf[cell[c->index()]][f->localIndex()];

                                // dof-id of other
                                BdryIdMapIterator j =
                                    bdry_id_map_.find(bc.getPeriodicPartner(bid));
                                ASSERT (j != bdry_id_map_.end());
                                const int c2   = j->second.first;
                                const int dof2 = cf[c2][j->second.second];

                                if (dof1 < dof2) {
                                    // Allocate storage for the actual
                                    // couplings.
                                    //
                                    const int ndof = cf.rowSize(c2);
                                    S_.incrementrowsize(dof1, ndof); // self->other
                                    for (int dof = 0; dof < ndof; ++dof) {
                                        int ii = cf[c2][dof];
                                        int pp = ppartner_dof_[ii];
                                        if ((pp != -1) && (pp != dof1) && (pp < ii)) {
                                            S_.incrementrowsize(pp, 1);
                                        }
                                        S_.incrementrowsize(ii, 1);  // other->self
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }



        // ----------------------------------------------------------------
        void setConnections(const BCInterface& bc)
        // ----------------------------------------------------------------
        {
            setGridConnections();
            setBCConnections(bc);

            S_.endindices();

            const int nc = pgrid_->numberOfCells();
            std::vector<Scalar>(nc).swap(flowSolution_.pressure_);
            std::vector<Scalar>(nc).swap(g_);
            std::vector<Scalar>(nc).swap(L_);

            matrix_structure_valid_ = true;
        }


        // ----------------------------------------------------------------
        void setGridConnections()
        // ----------------------------------------------------------------
        {
            const   SparseTable<int>& cf = flowSolution_.cellFaces_;
            typedef SparseTable<int>::row_type::const_iterator fi;

            // Compute actual connections (the non-zero structure).
            for (int c = 0; c < cf.size(); ++c) {
                fi fb = cf[c].begin(), fe = cf[c].end();

                for (fi i = fb; i != fe; ++i) {
                    for (fi j = fb; j != fe; ++j) {
                        S_.addindex(*i, *j);
                    }
                }
            }
        }


        // ----------------------------------------------------------------
        void setBCConnections(const BCInterface& bc)
        // ----------------------------------------------------------------
        {
            // Include system connections for periodic boundary
            // conditions, if any.  We make an arbitrary choice in
            // that the face/degree-of-freedom with the minimum
            // numerical id of the two periodic partners represents
            // the coupling.  Suppose <i_p> is this minimum dof-id.
            // We then need to introduce a *symmetric* coupling to
            // <i_p> to each of the dof's of the cell *NOT* containing
            // <i_p>.  This choice is implemented in the following
            // loop by introducing couplings only when dof1 (self) is
            // less than dof2 (other).
            //
            // See also: allocateBCConnections() and addCellContrib().
            //
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            if (!bdry_id_map_.empty()) {
                // At least one periodic BC.  Assign periodic
                // connections.
                for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
		    if ((*ppartition_)[c->index()] != my_partition_) continue;
                    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                        if (f->boundary()) {
                            const int bid = f->boundaryId();
                            if (bc.flowCond(bid).isPeriodic()) {
                                // dof-id of self
                                const int dof1 = cf[cell[c->index()]][f->localIndex()];

                                // dof-id of other
                                BdryIdMapIterator j =
                                    bdry_id_map_.find(bc.getPeriodicPartner(bid));
                                ASSERT (j != bdry_id_map_.end());
                                const int c2   = j->second.first;
                                const int dof2 = cf[c2][j->second.second];

                                if (dof1 < dof2) {
                                    // Assign actual couplings.
                                    //
                                    const int ndof = cf.rowSize(c2);
                                    for (int dof = 0; dof < ndof; ++dof) {
                                        int ii = cf[c2][dof];
                                        int pp = ppartner_dof_[ii];
                                        if ((pp != -1) && (pp != dof1) && (pp < ii)) {
                                            ii = pp;
                                        }
                                        S_.addindex(dof1, ii);  // self->other
                                        S_.addindex(ii, dof1);  // other->self
                                        S_.addindex(dof2, ii);
                                        S_.addindex(ii, dof2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }



        // ----------------------------------------------------------------
        void assembleDynamic(const ReservoirInterface&  r  ,
                             const std::vector<double>& sat,
                             const BCInterface&         bc ,
                             const std::vector<double>& src,
                             const typename GridInterface::CellIterator::Vector& grav)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            std::vector<double> mob(ReservoirInterface::NumberOfPhases);
            std::vector<double> rho(ReservoirInterface::NumberOfPhases);

            std::vector<Scalar> data_store(max_ncf_ * max_ncf_);
            std::vector<Scalar> e  (max_ncf_);
            std::vector<Scalar> rhs(max_ncf_);

            std::vector<FaceType> facetype(max_ncf_);
            std::vector<Scalar>   condval (max_ncf_);
            std::vector<int>      ppartner(max_ncf_);

            // Clear residual data
            S_   = 0.0;
            rhs_ = 0.0;

            std::fill(g_.begin(), g_.end(), Scalar(0.0));
            std::fill(L_.begin(), L_.end(), Scalar(0.0));
            std::fill(e .begin(), e .end(), Scalar(1.0));

            // We will have to regularize resulting system if there
            // are no prescribed pressures (i.e., Dirichlet BC's).
            do_regularization_ = true;

            typedef InnerProduct<CI,GridInterface::Dimension,true> IP;
            IP ip(max_ncf_);

            // Assemble dynamic contributions for each cell
            for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
                const int ci = c->index();
		if ((*ppartition_)[ci] != my_partition_) continue;
                const int c0 = cell[ci];            ASSERT (c0 < cf.size());
                const int nf = cf[c0].size();

                r.phaseMobility(ci, sat[ci], mob);
                r.phaseDensity (ci,          rho);

                const double totmob = std::accumulate   (mob.begin(), mob.end(), 0.0);
                const double omega  = std::inner_product(rho.begin(), rho.end(),
                                                         mob.begin(), 0.0) / totmob;

                SharedFortranMatrix    S  (nf, nf, &data_store[0]);
                ImmutableFortranMatrix one(nf, 1 , &e[0]);

                typename SparseTable<double>::mutable_row_type gterm = f_[c0];
                std::fill(gterm.begin(), gterm.end(), Scalar(0.0));
                ip.gravityTerm(c, grav, omega, gterm);

                setExternalContrib(c, c0, bc, src[ci], rhs,
                                   facetype, condval, ppartner);

                buildCellContrib(c0, totmob, one, S, rhs);

                addCellContrib(S, rhs, facetype, condval, ppartner, cf[c0]);
            }
        }



        // ----------------------------------------------------------------
        void solveLinearSystem(double residual_tolerance, int verbosity_level)
        // ----------------------------------------------------------------
        {
            // Adapted from DuMux...
            Scalar residTol = residual_tolerance;

            typedef BCRSMatrix <MatrixBlockType>        Matrix;
            typedef BlockVector<VectorBlockType>        Vector;
            typedef MatrixAdapter<Matrix,Vector,Vector> Adapter;

            // Regularize the matrix (only for pure Neumann problems...)
            if (do_regularization_) {
                S_[0][0] *= 2;
            }
            Adapter opS(S_);

            // Construct preconditioner.
            Dune::SeqILU0<Matrix,Vector,Vector> precond(S_, 1.0);

            // Construct solver for system of linear equations.
            Dune::BiCGSTABSolver<Vector> linsolve(opS, precond, residTol,
                                                  S_.N(), verbosity_level);

            Dune::InverseOperatorResult result;
            soln_ = 0.0;

            // Solve system of linear equations to recover
            // face/contact pressure values (soln_).
            linsolve.apply(soln_, rhs_, result);
        }


        /// @brief
        ///    Solve assemble system of linear equations @f$S\pi=r@f$
        ///    to recover interface/contact pressures, @f$\pi@f$.
        ///
        /// @param [in] residual_tolerance
        ///    Control parameter for iterative linear solver software.
        ///    The iteration process is terminated when the norm of
        ///    the linear system residual is less than @code
        ///    residual_tolerance @endcode times the initial residual.
        ///
        /// @param [in] linsolver_verbosity
        ///    Level of diagnostic output from linear system solver
        ///    output.  Semantics and possible values strongly tied to
        ///    the specifics of the ISTL PaAMG preconditioner and
        ///    CGSolver software.
        void solveLinearSystemAMG(double residual_tolerance,
                                  int verbosity_level)
        // ----------------------------------------------------------------
        {
            // Adapted from upscaling.cc by Arne Rekdal, 2009
            Scalar residTol = residual_tolerance;

            // Representation types for linear system.
            typedef BCRSMatrix <MatrixBlockType>        Matrix;
            typedef BlockVector<VectorBlockType>        Vector;
            // typedef MatrixAdapter<Matrix,Vector,Vector> Operator;
	    typedef OwnerOverlapCopyCommunication<LocalId,GlobalId> Communication;
	    typedef OverlappingSchwarzOperator<Matrix,Vector,Vector,Communication> Operator;  
	    typedef OverlappingSchwarzScalarProduct<Vector,Communication> ScalarProduct;

            // AMG specific types.
#define FIRST_DIAGONAL 1
#define SYMMETRIC 1

#if FIRST_DIAGONAL
            typedef Amg::FirstDiagonal CouplingMetric;
#else
            typedef Amg::RowSum        CouplingMetric;
#endif

#if SYMMETRIC
            typedef Amg::SymmetricCriterion<Matrix,CouplingMetric>   CriterionBase;
#else
            typedef Amg::UnSymmetricCriterion<Matrix,CouplingMetric> CriterionBase;
#endif

            typedef SeqILU0<Matrix, Vector, Vector>        Smoother;
	    typedef BlockPreconditioner<Vector, Vector, Communication, Smoother> ParSmoother;
            typedef Amg::CoarsenCriterion<CriterionBase> Criterion;
            typedef Amg::AMG<Operator, Vector, ParSmoother, Communication>   Precond;

            // Regularize the matrix (only for pure Neumann problems...)
            if (do_regularization_) {
                S_[0][0] *= 2;
            }
            Operator opS(S_, comm_);
	    ScalarProduct sp(comm_);

            // Construct preconditioner.
            double relax = 1;
            typename Precond::SmootherArgs smootherArgs;
            smootherArgs.relaxationFactor = relax;

            Criterion criterion;
            criterion.setDebugLevel(verbosity_level);
            Precond precond(opS, criterion, smootherArgs, 1, 2, 2, false, comm_);

            // Construct solver for system of linear equations.
            CGSolver<Vector> linsolve(opS, sp, precond, residTol, S_.N(), verbosity_level /*(rank == 0) ? verbosity_level : 0*/);

            InverseOperatorResult result;
            soln_ = 0.0;

            // Solve system of linear equations to recover
            // face/contact pressure values (soln_).
            linsolve.apply(soln_, rhs_, result);
        }



        /// @brief
        ///    Recover cell pressure and out-fluxes from
        ///    interface/contact pressure values.
        ///
        /// @param [in] r
        ///    The reservoir properties of each grid cell.
        ///
        /// @param [in] sat
        ///    Saturation of primary phase.  One scalar value for each
        ///    grid cell.  This parameter currently limits @code
        ///    IncompFlowSolverHybrid @endcode to two-phase flow
        ///    problems.
        void computePressureAndFluxes(const ReservoirInterface&  r  ,
                                      const std::vector<double>& sat)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            std::vector<Scalar>& p = flowSolution_.pressure_;
            SparseTable<Scalar>& v = flowSolution_.outflux_;

            std::vector<double> mob(ReservoirInterface::NumberOfPhases);
            std::vector<double> pi (max_ncf_);

            // Assemble dynamic contributions for each cell
            for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
                const int ci = c->index();
		if ((*ppartition_)[ci] != my_partition_) continue;
                const int c0 = cell[ci];
                const int nf = cf[c0].size();

                // Extract contact pressures for cell 'c'.
                for (int i = 0; i < nf; ++i) {
                    pi[i] = soln_[cf[c0][i]];
                }

                // Compute cell pressure in cell 'c'.
                p[c0] = (g_[c0] +
                         std::inner_product(F_[c0].begin(), F_[c0].end(),
                                            pi.begin(), 0.0)) / L_[c0];

                // Compute cell (out) fluxes for cell 'c'.
                // 1) Form system right hand side, r = f + Cp - D\pi
                std::transform(f_[c0].begin(), f_[c0].end(), pi.begin(),
                               pi    .begin(), //_1 + p[c0] - _2);
                               boost::bind(std::minus<Scalar>(),
                                           boost::bind(std::plus<Scalar>(),
                                                       _1,
                                                       p[c0]),
                                           _2));

                // 2) Solve system Bv = r
                r.phaseMobility(ci, sat[ci], mob);
                const double totmob = std::accumulate(mob.begin(), mob.end(), 0.0);

                ImmutableFortranMatrix Binv(nf, nf, &Binv_[c0][0]);
                vecMulAdd_N(totmob, Binv, &pi[0], Scalar(0.0), &v[c0][0]);
            }
        }




        // ----------------------------------------------------------------
        void setExternalContrib(const typename GridInterface::CellIterator c,
                                const int c0, const BCInterface& bc,
                                const double           src,
                                std::vector<Scalar>&   rhs,
                                std::vector<FaceType>& facetype,
                                std::vector<double>&   condval,
                                std::vector<int>&      ppartner)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator::FaceIterator FI;

            const SparseTable<int>& cf = flowSolution_.cellFaces_;

            std::fill(rhs     .begin(), rhs     .end(), Scalar(0.0));
            std::fill(facetype.begin(), facetype.end(), Internal   );
            std::fill(condval .begin(), condval .end(), Scalar(0.0));
            std::fill(ppartner.begin(), ppartner.end(), -1         );

            g_[c0] = src;

            int k = 0;
            for (FI f = c->facebegin(); f != c->faceend(); ++f, ++k) {
                if (f->boundary()) {
                    const int bid = f->boundaryId();
                    const FlowBC& bcond = bc.flowCond(bid);
                    if (bcond.isDirichlet()) {
                        facetype[k]        = Dirichlet;
                        condval[k]         = bcond.pressure();
                        do_regularization_ = false;
                    } else if (bcond.isPeriodic()) {
                        BdryIdMapIterator j =
                            bdry_id_map_.find(bc.getPeriodicPartner(bid));
                        ASSERT (j != bdry_id_map_.end());

                        facetype[k] = Periodic;
                        condval[k]  = bcond.pressureDifference();
                        ppartner[k] = cf[j->second.first][j->second.second];
                    } else {
                        ASSERT (bcond.isNeumann());
                        facetype[k] = Neumann;
                        rhs[k]      = bcond.outflux();
                    }
                }
            }
        }




        // ----------------------------------------------------------------
        void buildCellContrib(const int c, const Scalar totmob,
                              const ImmutableFortranMatrix& one,
                              SharedFortranMatrix& S, std::vector<Scalar>& rhs)
        // ----------------------------------------------------------------
        {
            std::transform(Binv_[c].begin(), Binv_[c].end(), S.data(),
                           boost::bind(std::multiplies<Scalar>(), _1, totmob));

            // Ft <- B^{-t} * ones([size(S,2),1])
            SharedFortranMatrix Ft(S.numRows(), 1, &F_[c][0]);
            matMulAdd_TN(Scalar(1.0), S, one, Scalar(0.0), Ft);

            L_[c]  = std::accumulate   (Ft.data(), Ft.data() + Ft.numRows(), 0.0);
            g_[c] -= std::inner_product(Ft.data(), Ft.data() + Ft.numRows(),
                                        f_[c].begin(), Scalar(0.0));

            // rhs <- B^{-1}*f - r (==B^{-1}f + E\pi - h)
            vecMulAdd_N(Scalar(1.0), S, &f_[c][0], -Scalar(1.0), &rhs[0]);

            // rhs <- rhs + g_[c]/L_[c]*F
            std::transform(rhs.begin(), rhs.end(), Ft.data(), rhs.begin(),
                           axpby<Scalar>(Scalar(1.0), Scalar(g_[c] / L_[c])));

            // S <- S - F'*F/L_c
            symmetricUpdate(-1.0/L_[c], Ft, 1.0, S);
        }



        /// @brief
        ///    Assemble local system contributions @f$S\pi = r@f$ of
        ///    single cell into global system.
        ///
        /// @param [in] S
        ///    Coefficient matrix contribution of single cell.
        ///
        /// @param [in] rhs
        ///    System right hand side contribution of single cell.
        ///
        /// @param [in] facetype
        ///    Indicator of which kind of equation to assemble for any
        ///    given interface/contact {Internal, Neumann, Dirichlet
        ///    or Periodic}.
        ///
        /// @param [in] condval
        ///    Numeric value of any external condition (e.g.,
        ///    prescribed face pressure for Dirichlet conditions,
        ///    prescribed fluxes for Neumann conditions, or prescribed
        ///    pressure drops for periodic conditions).
        ///
        /// @param [in] ppartner
        ///    A given DOF's periodic partner DOF.  Only examined in
        ///    the case of periodic boundary conditions.
        ///
        /// @param [in] l2g
        ///    Local to global face map for a single cell.
        template<class L2G>
        void addCellContrib(const SharedFortranMatrix&   S       ,
                            const std::vector<Scalar>&   rhs     ,
                            const std::vector<FaceType>& facetype,
                            const std::vector<Scalar>&   condval ,
                            const std::vector<int>&      ppartner,
                            const L2G&                   l2g)
        // ----------------------------------------------------------------
        {
            typedef typename L2G::const_iterator it;

            int r = 0;
            for (it i = l2g.begin(); i != l2g.end(); ++i, ++r) {
                // Indirection for periodic BC handling.
                int ii = *i;

                switch (facetype[r]) {
                case Dirichlet:
                    // Pressure is already known.  Assemble trivial
                    // equation of the form: a*x = a*p where 'p' is
                    // the known pressure value (i.e., condval[r]).
                    //
                    S_  [ii][ii] = S(r,r);
                    rhs_[ii]     = S(r,r) * condval[r];
                    break;
                case Periodic:
                    // Periodic boundary condition.  Contact pressures
                    // linked by relations of the form
                    //
                    //      a*(x0 - x1) = a*pd
                    //
                    // where 'pd' is the known pressure difference
                    // x0-x1 (condval[r]).  Preserve matrix symmetry
                    // by assembling both of the equations
                    //
                    //      a*(x0-x1) = a*  pd, and  (1)
                    //      a*(x1-x0) = a*(-pd)      (2)
                    //
                    ASSERT ((0 <= ppartner[r]) && (ppartner[r] < int(rhs_.size())));
                    ASSERT (ii != ppartner[r]);

                    {
                        const double a = S(r,r), b = a * condval[r];

                        // Equation (1)
                        S_  [         ii][         ii] += a;
                        S_  [         ii][ppartner[r]] -= a;
                        rhs_[         ii]              += b;

                        // Equation (2)
                        S_  [ppartner[r]][         ii] -= a;
                        S_  [ppartner[r]][ppartner[r]] += a;
                        rhs_[ppartner[r]]              -= b;
                    }

                    ii = std::min(ii, ppartner[r]);

                    // INTENTIONAL FALL-THROUGH!
                    // IOW: Don't insert <break;> here!
                    //
                default:
                    int c = 0;
                    for (it j = l2g.begin(); j != l2g.end(); ++j, ++c) {
                        // Indirection for periodic BC handling.
                        int jj = *j;

                        if (facetype[c] == Dirichlet) {
                            rhs_[ii] -= S(r,c) * condval[c];
                            continue;
                        }
                        if (facetype[c] == Periodic) {
                            ASSERT ((0 <= ppartner[c]) && (ppartner[c] < int(rhs_.size())));
                            ASSERT (jj != ppartner[c]);
                            if (ppartner[c] < jj) {
                                rhs_[ii] -= S(r,c) * condval[c];
                                jj = ppartner[c];
                            }
                        }
                        S_[ii][jj] += S(r,c);
                    }
                    break;
                }
                rhs_[ii] += rhs[r];
            }
        }
    };
} // namespace Dune

#endif // OPENRS_PARINCOMPFLOWSOLVERHYBRID_HEADER
