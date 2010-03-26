//===========================================================================
//
// File: known_answer_test.cpp
//
// Created: Thu Mar 25 13:57:12 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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



#include <config.h>

#include <algorithm>
#include <iostream>
#include <iomanip>

#include <boost/static_assert.hpp>

#include <dune/common/array.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/Units.hpp>

// #if HAVE_ALUGRID
// #include <dune/common/shared_ptr.hh>
// #include <dune/grid/io/file/gmshreader.hh>
// #include <dune/grid/alugrid.hh>
// #endif

#include <dune/solvers/common/SimulatorUtilities.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/common/EclipseGridParser.hpp>
#include <dune/grid/common/EclipseGridInspector.hpp>

#include <dune/solvers/common/fortran.hpp>
#include <dune/solvers/common/blas_lapack.hpp>
#include <dune/solvers/common/Matrix.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>

#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>
#include <dune/common/param/ParameterGroup.hpp>



// ------------ Specifying the solution ------------

typedef Dune::FieldVector<double, 3> Vec;

double u(const Vec& x)
{
    return std::sin(2*M_PI*x[0])*std::cos(2*M_PI*x[1])*x[2];
}
Vec Du(const Vec& x)
{
    Vec du;
    du[0] = 2*M_PI*std::cos(2*M_PI*x[0])*std::cos(2*M_PI*x[1])*x[2];
    du[1] = -2*M_PI*std::sin(2*M_PI*x[0])*std::sin(2*M_PI*x[1])*x[2];
    du[2] = 2*M_PI*std::sin(2*M_PI*x[0])*std::cos(2*M_PI*x[1]);
    return du;
}
double Lu(const Vec& x)
{
    return -2*2*M_PI*2*M_PI*std::sin(2*M_PI*x[0])*std::cos(2*M_PI*x[1])*x[2];
}

/*
double u(const Vec& x)
{
    return 0.5*x[0]*(1.0 - x[0]);
}
double Lu(const Vec& x)
{
    return -1.0;
}
*/
/*
double u(const Vec& x)
{
    return x[0]*x[1]*x[2];
}
Vec Du(const Vec& x)
{
    Vec du;
    du[0] = x[1]*x[2];
    du[1] = x[2]*x[0];
    du[2] = x[0]*x[1];
    return du;
}
double Lu(const Vec& x)
{
    return 0.0;
}
*/

/*
double u(const Vec& x)
{
    return x[0];
}
Vec Du(const Vec& x)
{
    Vec du;
    du[0] = 1.0;
    du[1] = 0.0;
    du[2] = 0.0;
    return du;
}
double Lu(const Vec& x)
{
    return 0.0;
}
*/


namespace Dune
{
    template <class BoundaryFunc>
    class FunctionBoundaryConditions : public PeriodicConditionHandler
    {
    public:
        FunctionBoundaryConditions(BoundaryFunc bfunc)
            : bfunc_(bfunc)
        {
        }

        template <class BoundaryFace>
        FlowBC flowCond(const BoundaryFace& bf) const
        {
            ASSERT(bf.boundary());
            return FlowBC(FlowBC::Dirichlet, bfunc_(bf.centroid()));
        }

    private:
        BoundaryFunc bfunc_;
    };

}

template<class GI>
void assign_src(const GI& g, std::vector<double>& src)
{
    typedef typename GI::CellIterator CI;
    int count = 0;
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
        src[count++] = -Lu(c->centroid())*c->volume();
    }
}

template<class GI, class BCS>
void assign_bc(const GI& g, BCS& bcs)
{
    typedef Dune::FlowBC BC;
    typedef typename GI::CellIterator CI;
    typedef typename CI::FaceIterator FI;
    int max_bid = 0;
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
        for (FI f = c->facebegin(); f != c->faceend(); ++f) {
            int bid = f->boundaryId();
            if (bid > max_bid) {
                max_bid = bid;
                bcs.resize(bid + 1);
            }
            bcs.flowCond(bid) = BC(BC::Dirichlet, u(f->centroid()));
        }
    }
}

template<class GI>
void compare_pressure(const GI& g, const std::vector<double>& p)
{
    typedef typename GI::CellIterator CI;
    int count = 0;
    double l1err = 0.0;
    double l2err = 0.0;
    double linferr = 0.0;
    double totv = 0.0;
    for (CI c = g.cellbegin(); c != g.cellend(); ++c, ++count) {
        Vec cen = c->centroid();
        double uval = u(cen);
        double diff = uval - p[count];
        double v = c->volume();
        l1err += std::fabs(diff*v);
        l2err += diff*diff*v;
        linferr = std::max(std::fabs(diff), linferr);
        totv += v;
        std::cout << cen[0] << ' ' << uval << ' ' << p[count] << std::endl;
    }
    l2err = std::sqrt(l2err);
    std::cout << "\n\n"
              << "\n     L1 error density: " << l1err/totv
              << "\n     L2 error density: " << l2err/totv
              << "\n     Linf error:       " << linferr << "\n\n\n";
}


template<class GI, class RI>
void test_flowsolver(const GI& g, const RI& r, double tol, int kind)
{
    typedef typename GI::CellIterator                   CI;
    typedef typename CI::FaceIterator                   FI;
    typedef double (*SolutionFuncPtr)(const Vec&);

    //typedef Dune::BasicBoundaryConditions<true, false>  FBC;
    typedef Dune::FunctionBoundaryConditions<SolutionFuncPtr> FBC;
    typedef Dune::IncompFlowSolverHybrid<GI, RI, FBC,
        Dune::MimeticIPEvaluator> FlowSolver;

    FlowSolver solver;

    // FBC flow_bc;
    // assign_bc(g, flow_bc);
    FBC flow_bc(&u);

    typename CI::Vector gravity(0.0);

    solver.init(g, r, gravity, flow_bc);

    std::vector<double> src(g.numberOfCells(), 0.0);
    assign_src(g, src);
    std::vector<double> sat(g.numberOfCells(), 0.0);


    solver.solve(r, sat, flow_bc, src, tol, 3, kind);

    typedef typename FlowSolver::SolutionType FlowSolution;
    FlowSolution soln = solver.getSolution();

    std::vector<double> cell_velocity;
    estimateCellVelocity(cell_velocity, g, soln);
    std::vector<double> cell_pressure;
    getCellPressure(cell_pressure, g, soln);

    compare_pressure(g, cell_pressure);

    Dune::VTKWriter<typename GI::GridType::LeafGridView> vtkwriter(g.grid().leafView());
    vtkwriter.addCellData(cell_velocity, "velocity");
    vtkwriter.addCellData(cell_pressure, "pressure");
    vtkwriter.write("testsolution-" + boost::lexical_cast<std::string>(0),
                    Dune::VTKOptions::ascii);
}



int main(int argc, char** argv)
{
    Dune::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Make a grid
    typedef Dune::CpGrid Grid;
    Grid grid;
    grid.init(param);
    grid.setUniqueBoundaryIds(true);

    // Make the grid interface
    Dune::GridInterfaceEuler<Grid> g(grid);

    // Reservoir properties.
    Dune::ReservoirPropertyCapillary<Grid::dimension> res_prop;
    res_prop.init(g.numberOfCells(), 1.0, 1.0);
    res_prop.setViscosities(1.0, 1.0);
    // res_prop.setDensities(1.0, 1.0);

    test_flowsolver(g, res_prop,
                    param.getDefault("tolerance", 1e-8),
                    param.getDefault("linear_solver_type", 1));
}
