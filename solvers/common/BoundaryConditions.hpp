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
#include <ostream>
#include <dune/common/ErrorMacros.hpp>


namespace Dune
{

    /// @brief A class for building boundary conditions in a uniform way.
    /// @todo Doc me!
    class BCBase
    {
    public:
	/// @brief Enum for the allowed boundary condition types.
	/// So far, we support Dirichlet, Neumann and Periodic conditions.
	/// In this class, these are just tags, it's up to the code using it
	/// to attach meaning to them.
        enum BCType { Dirichlet, Neumann, Periodic };

	/// @brief Write type and value to an output stream.
	/// @tparam traits character type.
	/// @tparam traits character traits.
	/// @param os output stream.
	template<typename charT, class traits>
	void write(std::basic_ostream<charT,traits>& os) const
	{
	    os << "Type: " << type_ << "   Value: " << value_;
	}

    protected: // methods
	/// @brief Default constructor, that makes a Neumann condition with value 0.0.
        BCBase()
            : type_(Neumann), value_(0.0)
        {
        }
	/// @brief Constructor taking a type and value.
	/// @param type the condition type.
	/// @param value the condition value.
        BCBase(BCType type, double value)
            : type_(type), value_(value)
        {
        }
	/// @brief Type query.
	/// @return true if the type is Dirichlet.
        bool isDirichlet() const
        {
            return type_ == Dirichlet;
        }
	/// @brief Type query.
	/// @return true if the type is Neumann.
        bool isNeumann() const
        {
            return type_ == Neumann;
        }
	/// @brief Type query.
	/// @return true if the type is Periodic.
        bool isPeriodic() const
        {
            return type_ == Periodic;
        }

    protected: // data
	BCType type_;
	double value_;
    };

    template<typename charT, class traits>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const BCBase& bc)
    {
        bc.write(os);
        return os;
    }


    /// @brief A class for representing a flow boundary condition.
    /// @todo Doc me!
    class FlowBC : public BCBase
    {
    public:
	/// @brief Default constructor, that makes a Neumann condition with value 0.0.
        FlowBC()
            : BCBase()
        {
        }
	/// @brief Constructor taking a type and value.
	/// @param type the condition type.
	/// @param value the condition value.
        FlowBC(BCType type, double value)
            : BCBase(type, value)
        {
        }

	using BCBase::isDirichlet;
	using BCBase::isNeumann;
	using BCBase::isPeriodic;

	/// @brief Query a Dirichlet condition.
	/// @return the pressure condition value
        double pressure() const
        {
            ASSERT(type_ == Dirichlet);
            return value_;
        }
	/// @brief Query a Neumann condition.
	/// @return the outwards flux condition value.
        double outflux() const
        {
            ASSERT(type_ == Neumann);
            return value_;
        }
	/// @brief Query a Periodic condition.
	/// @return the pressure difference condition value.
        double pressureDifference() const
        {
            ASSERT(type_ == Periodic);
            return value_;
        }
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

        void resize(int new_size, FlowBC fbc)
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

        const FlowBC& operator[](int boundary_id) const
        {
            ASSERT(boundary_id >= 0 && boundary_id < int(fbcs_.size()));
            return fbcs_[boundary_id];
        }

        FlowBC& operator[](int boundary_id)
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

        template<typename charT, class traits>
        void write(std::basic_ostream<charT,traits>& os) const
        {
            for (int i = 0;  i < int(fbcs_.size()); ++i) {
                os << fbcs_[i] << "   " << periodic_partner_bid_[i] << '\n';
            }
            os << std::endl;
        }
    private:
        std::vector<FlowBC> fbcs_;
        std::vector<int> periodic_partner_bid_;
    };

    template<typename charT, class traits>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const FlowBoundaryConditions& fbcs)
    {
        fbcs.write(os);
        return os;
    }





    class SaturationBoundaryCondition
    {
    public:
        enum BCType { Dirichlet, Periodic };
        SaturationBoundaryCondition()
            : type_(Dirichlet), value_(1.0)
        {
        }
        SaturationBoundaryCondition(BCType type, double value)
            : type_(type), value_(value)
        {
        }
	/// @brief
	/// @todo Doc me!
	/// @return
        bool isDirichlet() const
        {
            return type_ == Dirichlet;
        }
	/// @brief
	/// @todo Doc me!
	/// @return
        bool isPeriodic() const
        {
            return type_ == Periodic;
        }
	/// @brief
	/// @todo Doc me!
	/// @return
        double saturation() const
        {
            ASSERT(type_ == Dirichlet);
            return value_;
        }
	/// @brief
	/// @todo Doc me!
	/// @return
        double saturationDifference() const
        {
            ASSERT(type_ == Periodic);
            return value_;
        }
	/// @brief
	/// @todo Doc me!
	/// @tparam
	/// @param
        template<typename charT, class traits>
        void write(std::basic_ostream<charT,traits>& os) const
        {
            os << "Type: " << type_ << "   Value: " << value_;
        }
    private:
        BCType type_;
        double value_;
    };

    template<typename charT, class traits>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const SaturationBoundaryCondition& sbc)
    {
        sbc.write(os);
        return os;
    }



    class SaturationBoundaryConditions
    {
    public:
        SaturationBoundaryConditions()
        {
        }

        SaturationBoundaryConditions(int num_different_boundary_ids)
            : sbcs_(num_different_boundary_ids),
              periodic_partner_bid_(num_different_boundary_ids, 0)
        {
        }

        void resize(int new_size)
        {
            sbcs_.resize(new_size);
            periodic_partner_bid_.resize(new_size, 0);
        }

        void resize(int new_size, SaturationBoundaryCondition sbc)
        {
            sbcs_.resize(new_size, sbc);
            periodic_partner_bid_.resize(new_size, 0);
        }

        bool empty() const
        {
            return sbcs_.empty();
        }

        void clear()
        {
            sbcs_.clear();
            periodic_partner_bid_.clear();
        }

        const SaturationBoundaryCondition& operator[](int boundary_id) const
        {
            ASSERT(boundary_id >= 0 && boundary_id < int(sbcs_.size()));
            return sbcs_[boundary_id];
        }

        SaturationBoundaryCondition& operator[](int boundary_id)
        {
            ASSERT(boundary_id >= 0 && boundary_id < int(sbcs_.size()));
            return sbcs_[boundary_id];
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

        template<typename charT, class traits>
        void write(std::basic_ostream<charT,traits>& os) const
        {
            for (int i = 0;  i < int(sbcs_.size()); ++i) {
                os << sbcs_[i] << "   " << periodic_partner_bid_[i] << '\n';
            }
            os << std::endl;
        }
    private:
        std::vector<SaturationBoundaryCondition> sbcs_;
        std::vector<int> periodic_partner_bid_;
    };

    template<typename charT, class traits>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const SaturationBoundaryConditions& sbcs)
    {
        sbcs.write(os);
        return os;
    }


} // namespace Dune


#endif // OPENRS_BOUNDARYCONDITIONS_HEADER
