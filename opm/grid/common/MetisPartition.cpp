/*
  Copyright 2024 OPM-OP AS

  This file is part of The Open Porous Media project  (OPM).

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_MPI // no code in this file without MPI, then skip includes.
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/MetisPartition.hpp>
#include <opm/grid/utility/OpmWellType.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <algorithm>
#include <type_traits>
#endif

#if defined(HAVE_METIS) && HAVE_MPI
namespace Dune
{
namespace cpgrid
{

// We want to use METIS, but if METIS is installed together with Scotch, then the following options are not available.
void setMetisOptions(const std::map<std::string, std::string>& optionsMap, idx_t* options) {

    // Initialize all options to default values
    METIS_SetDefaultOptions(options);

    // A map to translate string keys to METIS option indices
    // This is the list of options available for METIS Version 5.1.0 - possibly more can be added in the future
    std::map<std::string, int> metisOptionKeys = {
        // These options are only valid for the METIS_PartGraphKway method
        {"METIS_OPTION_OBJTYPE", METIS_OPTION_OBJTYPE},
        {"METIS_OPTION_MINCONN", METIS_OPTION_MINCONN},
        {"METIS_OPTION_CONTIG", METIS_OPTION_CONTIG},
        // The options are vaild for the METIS_PartGraphKway method and the METIS_PartGraphRecursive method
        {"METIS_OPTION_CTYPE", METIS_OPTION_CTYPE},
        {"METIS_OPTION_IPTYPE", METIS_OPTION_IPTYPE},
        {"METIS_OPTION_RTYPE", METIS_OPTION_RTYPE},
        {"METIS_OPTION_NO2HOP", METIS_OPTION_NO2HOP},
        {"METIS_OPTION_NCUTS", METIS_OPTION_NCUTS},
        {"METIS_OPTION_NITER", METIS_OPTION_NITER},
        {"METIS_OPTION_SEED", METIS_OPTION_SEED},
        {"METIS_OPTION_UFACTOR", METIS_OPTION_UFACTOR},
        {"METIS_OPTION_NUMBERING", METIS_OPTION_NUMBERING},
        {"METIS_OPTION_DBGLVL", METIS_OPTION_DBGLVL}
    };


    // Iterate over the input map and set the options accordingly
    for (const auto& pair : optionsMap) {
        const std::string& key = pair.first;
        const std::string& value = pair.second;

        if (metisOptionKeys.find(key) != metisOptionKeys.end()) {
            idx_t optionIndex = metisOptionKeys[key];
            options[optionIndex] = std::stoi(value); // Convert the value to integer and set the option
            Opm::OpmLog::info("Set metis option" + key + " to " + value + ".");
        } else {
            std::cerr << "Unknown METIS option: " << key << std::endl;
        }
    }
}


std::tuple<std::vector<int>,
           std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<int, int, char>>,
           std::vector<std::tuple<int, int, char, int>>,
           WellConnections>
metisSerialGraphPartitionGridOnRoot(const CpGrid& cpgrid,
                                    const std::vector<OpmWellType> * wells,
                                    const double* transmissibilities,
                                    const Communication<MPI_Comm>& cc,
                                    EdgeWeightMethod edgeWeightsMethod,
                                    int root,
                                    const real_t imbalanceTol,
                                    bool allowDistributedWells,
                                    [[maybe_unused]] const std::map<std::string,std::string>& params)
{
#if defined(IDXTYPEWIDTH) && IDXTYPEWIDTH != 64
    if (edgeWeightsMethod == Dune::EdgeWeightMethod::defaultTransEdgeWgt )
        OPM_THROW(std::runtime_error, "The selected partition method is METIS with default edge weights (i.e. transmissibilities).\
            This combination works only for METIS with 64-bit integers, but the installed version of METIS does not use 64-bit integers.\
            Either reinstall METIS with 64-bit integers or choose another edge weight method!");
#endif

    std::shared_ptr<CombinedGridWellGraph> gridAndWells;
    if( wells )
    {
        gridAndWells.reset(new CombinedGridWellGraph(cpgrid,
                                                       wells,
                                                       transmissibilities,
                                                       false,
                                                       edgeWeightsMethod));
    }

    std::vector<int> partitionVector;
    int rc = METIS_OK;

    cc.barrier();
    //Metis is a serial graph partitioner, we do everything only on root
    if (cc.rank() == root) {

        //////// First, we define all variables that *do not depend* on whether there are wells or not

        // The number of vertices, every cell is a vertex in the graph.
        idx_t n = cpgrid.numCells();
        
        // This is a vector of size n that upon successful completion stores the partition vector of the graph.
        // The numbering of this vector starts from either 0 or 1, depending on the value of options[METIS_OPTION_NUMBERING].
        idx_t* gpart = new idx_t[n];

        // Upon successful completion, this variable stores the edge-cut or the total communication volume of
        // the partitioning solution. The value returned depends on the partitioning’s objective function.
        idx_t objval = 0;

        // The number of partitions to split the graph into, we want to distribtue over all processes, so cc.size()
        idx_t nparts = cc.size(); 

        auto& globalIdSet         =  cpgrid.globalIdSet();
        auto& localIdSet          =  cpgrid.localIdSet();

        idx_t* gids = new idx_t[n];
        idx_t* lids = new idx_t[n];

        int idx = 0;
        for (auto cell = cpgrid.leafbegin<0>(), cellEnd = cpgrid.leafend<0>(); cell != cellEnd; ++cell)
        {
            gids[idx]   = globalIdSet.id(*cell);
            lids[idx++] = localIdSet.id(*cell);
        }

        //The number of balancing constraints, should be at least 1.
        idx_t ncon = 1;

        // The adjacency structure of a graph with n vertices and m edges is represented using two arrays xadj and adjncy.
        // An array of size n+1 that specifies the adjacency structure of the graph. The adjacency list of vertex i is stored in adjncy[xadj[i]] to adjncy[xadj[i+1]-1].
        idx_t* xadj = new idx_t[n+1];
        xadj[0] = 0;

        // This is an array of size ncon (in our case, of size 1) that specifies the allowed load imbalance tolerance for each constraint.
        // For the ith partition and jth constraint the allowed weight is the ubvec[j]*tpwgts[i*ncon+j] fraction
        // of the jth’s constraint total weight. The load imbalances must be greater than 1.0.
        // A NULL value can be passed indicating that the load imbalance tolerance for each constraint should
        // be 1.001 (for ncon=1) or 1.01 (for ncon>1).
        // NOTE: Scotch and METIS interpret this parameter differently
        real_t ubvec = imbalanceTol;
        
        // This is the array of options as described in Section 5.4.
        // The METIS options are not available if METIS is installed together with Scotch.
        idx_t* options = new idx_t[METIS_NOPTIONS];
        Dune::cpgrid::setMetisOptions(params, options);

        //////// Now, we define all variables that *do depend* on whether there are wells or not

        if( wells )
        {            
            for (int i = 0; i < n;  i++) {
                xadj[i+1] = xadj[i] + Dune::cpgrid::getNumberOfEdgesForSpecificCellForGridWithWells(*gridAndWells, lids[i]);
            }
        }
        else
        {
            for (int i = 0; i < n;  i++) {
                xadj[i+1] = xadj[i] + Dune::cpgrid::getNumberOfEdgesForSpecificCell(cpgrid, lids[i]);
            }
        }

        // The number of edges depends on whether there are wells or not, twoM = 2*m, where m = number of edges.
        idx_t twoM = xadj[n];

        // An array that contains the adjacency list of the graph.
        // The xadj array is of size n + 1 whereas the adjncy array is of size 2m (because for each edge between vertices v and u we actually store both (v, u) and (u, v)).
        // The adjacency list of vertex i is stored in array adjncy starting at index xadj[i] and ending at (but not
        // including) index xadj[i + 1] (i.e., adjncy[xadj[i]] through and including adjncy[xadj[i + 1]-1])
        // So: xadj contains the indices where we start for the respective component
        idx_t* adjncy = new idx_t[twoM]; 
        
        // An array that contains the weights of the edges. If all edges have the same weight, this can be set to NULL.
        // The weights of the edges (if any) are stored in an additional array called adjwgt. This array contains 2m elements, and the weight of edge adjncy[j] is stored at location adjwgt[j]
        idx_t* adjwgt = new idx_t[twoM];

        if( wells )
        {
            int neighborCounter = 0;
            for( int cell = 0; cell < n;  cell++ )
            {
                fillNBORGIDAndWeightsForSpecificCellAndIncrementNeighborCounterForGridWithWells(*gridAndWells, lids[cell], gids, neighborCounter, adjncy, adjwgt);
            }
        }
        else
        {
            int neighborCounter = 0;
            for( int cell = 0; cell < n;  cell++ )
            {
                fillNBORGIDForSpecificCellAndIncrementNeighborCounter(cpgrid, lids[cell], gids, neighborCounter, adjncy);
            }
        }

        // Decide which partition method to use, both methods create k partitions, where
        // METIS_PartGraphRecursive uses multilevel recursive bisection and
        // METIS_PartGraphKway uses multilevel k-way partition.
        // The advice is: Use METIS_PartGraphRecursive if k is small and if k is a power of two
        // (n & (n - 1) == 0) is true if n > 0 and n is a power of two, this is an efficient bitwise check.
        if (nparts < 65 && ((nparts & (nparts - 1)) == 0)) {
            Opm::OpmLog::info("Partitioning grid using METIS_PartGraphRecursive.");
            rc = METIS_PartGraphRecursive(&n, 
                                          &ncon,
                                          xadj,
                                          adjncy,
                                          nullptr, // vwgt
                                          nullptr, // vsize,
                                          wells ? adjwgt : nullptr,
                                          &nparts,
                                          nullptr, // tpwgts,
                                          &ubvec,
                                          options,
                                          &objval,
                                          gpart);
        } else {
            Opm::OpmLog::info("Partitioning grid using METIS_PartGraphKway.");
            rc = METIS_PartGraphKway(&n, 
                                     &ncon,
                                     xadj,
                                     adjncy,
                                     nullptr, // vwgt
                                     nullptr, // vsize,
                                     wells ? adjwgt : nullptr,
                                     &nparts,
                                     nullptr, // tpwgts,
                                     &ubvec,
                                     options,
                                     &objval,
                                     gpart);
        }

        partitionVector.assign(gpart, gpart + n);
        
        delete[] gids;
        delete[] lids;
        delete[] xadj;
        delete[] adjncy;
        delete[] adjwgt;
        delete[] options;
        delete[] gpart;
    }

    //Broadcast the return value to all processes
    cc.broadcast(&rc, 1, root);
    if (rc == METIS_OK) {
        // Function returned normally :)
    } else if (rc == METIS_ERROR_INPUT) {
        OPM_THROW(std::runtime_error, "METIS Input Error!");
    } else if (rc == METIS_ERROR_MEMORY) {
        OPM_THROW(std::runtime_error, "METIS could not allocate the required memory!");
    } else if (rc == METIS_ERROR) {
        OPM_THROW(std::runtime_error, "Some other type of METIS error!");
    } else {   
        OPM_THROW(std::runtime_error, "Some other type of general error!");
    }
    
    return cpgrid::createListsFromParts(cpgrid, wells, transmissibilities, partitionVector, allowDistributedWells, gridAndWells);
}

} // namespace cpgrid
} // namespace Dune
#endif // HAVE_METIS && HAVE_MPI
