/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2014 Statoil ASA.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulartion-Software & Services
  Copyright 2015       NTNU

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

#ifndef OPM_GRID_ENUMS_HPP
#define OPM_GRID_ENUMS_HPP

namespace Dune {
    /// \brief enum for choosing Methods for weighting graph-edges correspoding to cell interfaces in Zoltan's or Metis' graph partitioner.
    ////
    /// uniform methods means all edges have weight 1. defaultTrans uses transmissibility as weights.
    /// logTrans uses the logarithm of the transmissibility.
    /// The uniform and logTrans edge-weighting methods produce partitioning results with lower edge-cut,
    /// fewer overlap/ghost cells and less communication overhead than when using defaultTrans. However, the impact
    /// on parallel linear solver performance is negative.
    enum EdgeWeightMethod {
        /// \brief All edge have a uniform weight of 1
        uniformEdgeWgt=0,
        /// \brief Use the transmissibilities as edge weights
        defaultTransEdgeWgt=1,
        /// \brief Use the log of the transmissibilities as edge weights
        logTransEdgeWgt=2
    };

    /// \brief enum for choosing methods for partitioning a graph.
    enum PartitionMethod {
        /// \brief Use simple approach based on rectangular partitioning the underlying cartesian grid.
        simple=0,
        /// \brief Use Zoltan for partitioning
        zoltan=1,
        /// \brief Use METIS for partitioning
        metis=2,
        /// \brief use Zoltan on GraphOfGrid for partitioning
        zoltanGoG=3
    };
}

#endif
