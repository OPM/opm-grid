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
#include <boost/mpl/if.hpp>
#include <dune/common/ErrorMacros.hpp>

namespace Dune
{

    /// @brief A class for building boundary conditions in a uniform way.
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

    /// @brief Stream insertion for BCBase.
    template<typename charT, class traits>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const BCBase& bc)
    {
        bc.write(os);
        return os;
    }




    /// @brief A class for representing a flow boundary condition.
    class FlowBC : public BCBase
    {
    public:
	/// @brief Default constructor, that makes a noflow condition (Neumann, value 0.0).
        FlowBC()
            : BCBase(Neumann, 0.0)
        {
        }
	/// @brief Constructor taking a type and value.
	/// @param type the condition type.
	/// @param value the condition value.
        FlowBC(BCType type, double value)
            : BCBase(type, value)
        {
	    ASSERT(isNeumann() || isDirichlet() || isPeriodic());
        }

	/// @brief Forwarding the relevant type queries.
	using BCBase::isDirichlet;
	using BCBase::isNeumann;
	using BCBase::isPeriodic;

	/// @brief Query a Dirichlet condition.
	/// @return the pressure condition value
        double pressure() const
        {
            ASSERT (isDirichlet());
            return value_;
        }
	/// @brief Query a Neumann condition.
	/// @return the outwards flux condition value.
        double outflux() const
        {
            ASSERT (isNeumann());
            return value_;
        }
	/// @brief Query a Periodic condition.
	/// @return the pressure difference condition value.
        double pressureDifference() const
        {
            ASSERT (isPeriodic());
            return value_;
        }
    };



    /// @brief A class for representing a saturation boundary condition.
    class SatBC : public BCBase
    {
    public:
	/// @brief Default constructor, that makes a Dirichlet condition with value 1.0.
        SatBC()
            : BCBase(Dirichlet, 1.0)
        {
        }
	/// @brief Constructor taking a type and value.
	/// @param type the condition type.
	/// @param value the condition value.
        SatBC(BCType type, double value)
            : BCBase(type, value)
        {
	    ASSERT(isDirichlet() || isPeriodic());
        }
	/// @brief Forwarding the relevant type queries.
	using BCBase::isDirichlet;
	using BCBase::isPeriodic;

	/// @brief Query a Dirichlet condition.
	/// @return the boundary saturation value
        double saturation() const
        {
            ASSERT (isDirichlet());
            return value_;
        }
	/// @brief Query a Periodic condition.
	/// @return the saturation difference value.
        double saturationDifference() const
        {
            ASSERT (isPeriodic());
            return value_;
        }
    };


    /// @brief A class for representing a capillary pressure boundary condition.
    class PcapBC : public BCBase
    {
    public:
	/// @brief Default constructor, that makes a Dirichlet condition with value 0.0.
        PcapBC()
            : BCBase(Dirichlet, 0.0)
        {
        }
	/// @brief Constructor taking a type and value.
	/// @param type the condition type.
	/// @param value the condition value.
        PcapBC(BCType type, double value)
            : BCBase(type, value)
        {
	    ASSERT(isDirichlet() || isPeriodic());
        }
	/// @brief Forwarding the relevant type queries.
	using BCBase::isDirichlet;
	using BCBase::isPeriodic;

	/// @brief Query a Dirichlet condition.
	/// @return the boundary saturation value
        double capPressure() const
        {
            ASSERT (isDirichlet());
            return value_;
        }
	/// @brief Query a Periodic condition.
	/// @return the saturation difference value.
        double capPressureDifference() const
        {
            ASSERT (isPeriodic());
            return value_;
        }
    };




    class PeriodicConditionHandler
    {
    public:
        PeriodicConditionHandler()
        {
        }

        PeriodicConditionHandler(int num_different_boundary_ids)
	    : periodic_partner_bid_(num_different_boundary_ids, 0)
        {
        }

        void resize(int new_size)
        {
            periodic_partner_bid_.resize(new_size, 0);
        }

        bool empty() const
        {
            return periodic_partner_bid_.empty();
        }

        void clear()
        {
            periodic_partner_bid_.clear();
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
            for (int i = 0;  i < int(periodic_partner_bid_.size()); ++i) {
                os << "Partner of bid " << i << " is " << periodic_partner_bid_[i] << '\n';
            }
            os << std::endl;
        }
    private:
        std::vector<int> periodic_partner_bid_;
    };

    template<typename charT, class traits>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const PeriodicConditionHandler& pch)
    {
        pch.write(os);
        return os;
    }


    template <typename T>
    class DummyVec
    {
    public:
	DummyVec() {}
	DummyVec(int) {}
	void resize(int) {}
	void clear() {}
    };

    template <bool FC = false, bool SC = false, bool PC = false >
    class BoundaryConditions : public PeriodicConditionHandler,
			       private boost::mpl::if_c<FC, std::vector<FlowBC>, DummyVec<FlowBC> >::type,
			       private boost::mpl::if_c<SC, std::vector<SatBC>,  DummyVec<SatBC>  >::type,
			       private boost::mpl::if_c<PC, std::vector<PcapBC>, DummyVec<PcapBC> >::type
    {
    public:
	typedef typename boost::mpl::if_c<FC, std::vector<FlowBC>, DummyVec<FlowBC> >::type FlowConds;
	typedef typename boost::mpl::if_c<SC,  std::vector<SatBC>, DummyVec<SatBC>  >::type SatConds;
	typedef typename boost::mpl::if_c<PC, std::vector<PcapBC>, DummyVec<PcapBC> >::type PcapConds;
	const static bool HasFlowConds = FC;
	const static bool HasSatConds = SC;
	const static bool HasPcapConds = PC;

	FlowBC& flowCond(int i)
	{
	    return FlowConds::operator[](i);
	}
	const FlowBC& flowCond(int i) const
	{
	    return FlowConds::operator[](i);
	}
	SatBC& satCond(int i)
	{
	    return SatConds::operator[](i);
	}
	const SatBC& satCond(int i) const
	{
	    return SatConds::operator[](i);
	}
	PcapBC& pcapCond(int i)
	{
	    return PcapConds::operator[](i);
	}
	const PcapBC& pcapCond(int i) const
	{
	    return PcapConds::operator[](i);
	}

        BoundaryConditions()
	{
        }

        BoundaryConditions(int num_different_boundary_ids)
	    : PeriodicConditionHandler(num_different_boundary_ids),
	      FlowConds(num_different_boundary_ids),
	      SatConds(num_different_boundary_ids),
	      PcapConds(num_different_boundary_ids)
	{
        }

        void resize(int new_size)
        {
	    PeriodicConditionHandler::resize(new_size);
	    FlowConds::resize(new_size);
	    SatConds::resize(new_size);
	    PcapConds::resize(new_size);
        }

        bool empty() const
        {
            return PeriodicConditionHandler::empty();
        }

        void clear()
        {
	    PeriodicConditionHandler::clear();
	    FlowConds::clear();
	    SatConds::clear();
	    PcapConds::clear();
        }

        template<typename charT, class traits>
        void write(std::basic_ostream<charT,traits>& os) const
        {
	    PeriodicConditionHandler::write(os);
	}
    };

    template<typename charT, class traits, bool F, bool S, bool P>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const BoundaryConditions<F,S,P>& bcs)
    {
        bcs.write(os);
        return os;
    }


} // namespace Dune


#endif // OPENRS_BOUNDARYCONDITIONS_HEADER
