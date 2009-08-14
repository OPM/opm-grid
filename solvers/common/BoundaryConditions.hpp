//===========================================================================
//
// File: BoundaryConditions.hpp
//
// Created: Mon Jun 29 15:19:42 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPENRS_BOUNDARYCONDITIONS_HEADER
#define OPENRS_BOUNDARYCONDITIONS_HEADER


#include <vector>

#include <dune/common/ErrorMacros.hpp>


namespace Dune
{

    class FlowBoundaryCondition
    {
    public:
        enum BCType { Dirichlet, Neumann, Periodic };
        FlowBoundaryCondition()
            : type_(Neumann), value_(0.0)
        {
        }
        FlowBoundaryCondition(BCType type, double value)
            : type_(type), value_(value)
        {
        }
        bool isDirichlet() const
        {
            return type_ == Dirichlet;
        }
        bool isNeumann() const
        {
            return type_ == Neumann;
        }
        bool isPeriodic() const
        {
            return type_ == Periodic;
        }
        double pressure() const
        {
            ASSERT(type_ == Dirichlet);
            return value_;
        }
        double outflux() const
        {
            ASSERT(type_ == Neumann);
            return value_;
        }
        double pressureDifference() const
        {
            ASSERT(type_ == Periodic);
            return value_;
        }
    private:
        BCType type_;
        double value_;
    };

    class FlowBoundaryConditions
    {
    public:
        FlowBoundaryConditions()
        {
        }

        FlowBoundaryConditions(int num_different_boundary_ids)
            : fbcs_(num_different_boundary_ids),
              periodic_partner_bid_(num_different_boundary_ids, 0)
        {
        }

        void resize(int new_size)
        {
            fbcs_.resize(new_size);
            periodic_partner_bid_.resize(new_size, 0);
        }

        void resize(int new_size, FlowBoundaryCondition fbc)
        {
            fbcs_.resize(new_size, fbc);
            periodic_partner_bid_.resize(new_size, 0);
        }

        bool empty() const
        {
            return fbcs_.empty();
        }

        void clear()
        {
            fbcs_.clear();
	    periodic_partner_bid_.clear();
        }

        const FlowBoundaryCondition& operator[](int boundary_id) const
        {
            ASSERT(boundary_id >= 0 && boundary_id < int(fbcs_.size()));
            return fbcs_[boundary_id];
        }

        FlowBoundaryCondition& operator[](int boundary_id)
        {
            ASSERT(boundary_id >= 0 && boundary_id < int(fbcs_.size()));
            return fbcs_[boundary_id];
        }

        void setPeriodicPartners(int boundary_id_1, int boundary_id_2)
        {
            ASSERT(boundary_id_1 >= 0 && boundary_id_1 < int(periodic_partner_bid_.size()));
            ASSERT(boundary_id_2 >= 0 && boundary_id_2 < int(periodic_partner_bid_.size()));
            ASSERT(periodic_partner_bid_[boundary_id_1] == 0);
            ASSERT(periodic_partner_bid_[boundary_id_2] == 0);
            periodic_partner_bid_[boundary_id_1] = boundary_id_2;
            periodic_partner_bid_[boundary_id_2] = boundary_id_1;
        }

        int getPeriodicPartner(int boundary_id) const
        {
            ASSERT(boundary_id >= 0 && boundary_id < int(periodic_partner_bid_.size()));
            return periodic_partner_bid_[boundary_id];
        }

    private:
        std::vector<FlowBoundaryCondition> fbcs_;
        std::vector<int> periodic_partner_bid_;
    };


    class SaturationBoundaryCondition
    {
    public:
        enum BCType { Dirichlet };
        SaturationBoundaryCondition()
            : type_(Dirichlet), value_(1.0)
        {
        }
        SaturationBoundaryCondition(BCType type, double value)
            : type_(type), value_(value)
        {
        }
        bool isDirichlet() const
        {
            return type_ == Dirichlet;
        }
        double saturation() const
        {
            ASSERT(type_ == Dirichlet);
            return value_;
        }
    private:
        BCType type_;
        double value_;
    };

    typedef std::vector<SaturationBoundaryCondition> SaturationBoundaryConditions;


} // namespace Dune


#endif // OPENRS_BOUNDARYCONDITIONS_HEADER
