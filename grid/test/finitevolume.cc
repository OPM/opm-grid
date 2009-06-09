#include "config.h"               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/common/mpihelper.hh> // include mpi helper class

// checks for defined gridtype and inlcudes appropriate dgfparser implementation
//#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include "vtkout.hh"
// #include"unitcube.hh"
// #include"transportproblem2.hh"
#include "initialize.hh"
#include "evolve.hh"

#include "../CpGrid.hpp"

typedef Dune::CpGrid GridType;

//===============================================================
// the time loop function working for all types of grids
//===============================================================

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
    bool contains(Dune::GeometryType gt)
    {
        return gt.dim() == dim;
    }
};

template<class G>
void timeloop(const G& grid, double tend)
{
    // make a mapper for codim 0 entities in the leaf grid
    Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout>
    mapper(grid);

    // allocate a vector for the concentration
    std::vector<double> c(mapper.size());

    // initialize concentration with initial values
    initialize(grid,mapper,c);                           /*@\label{fvc:init}@*/

    vtkout(grid,c,"concentration",0,0.0);

    // now do the time steps
    double t=0,dt;
    int k=0;
    const double saveInterval = 0.1;
    double saveStep = 0.1;
    int counter = 1;

    while (t<tend)                                       /*@\label{fvc:loop0}@*/
    {
        // augment time step counter
        ++k;

        // apply finite volume scheme
        evolve(grid,mapper,c,t,dt);

        // augment time
        t += dt;

        // check if data should be written
        if (t >= saveStep)
        {
            // write data
            vtkout(grid,c,"concentration",counter,t);

            // increase counter and saveStep for next interval
            saveStep += saveInterval;
            ++counter;
        }

        // print info about time, timestep size and counter
        std::cout << "s=" << grid.size(0)
        << " k=" << k << " t=" << t << " dt=" << dt << std::endl;
    }                                              /*@\label{fvc:loop1}@*/

    // output results
    vtkout(grid,c,"concentration",counter,tend);     /*@\label{fvc:file}@*/
}


//===============================================================
// The main function creates objects and does the time loop
//===============================================================

int main(int argc , char ** argv)
{
    // initialize MPI, finalize is done automatically on exit
    Dune::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // start try/catch block to get error messages from dune
    try {
        using namespace Dune;

        GridType grid;
	grid.init(param);

        // do time loop until end time 0.5
        timeloop(grid, 0.5);
    } catch (std::exception & e) {
        std::cout << "STL ERROR: " << e.what() << std::endl;
        return 1;
    } catch (Dune::Exception & e) {
        std::cout << "DUNE ERROR: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cout << "Unknown ERROR" << std::endl;
        return 1;
    }
    // done
    return 0;
}
