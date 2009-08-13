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
	enum BCType { Dirichlet, Neumann };
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
    private:
	BCType type_;
	double value_;
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


    typedef std::vector<FlowBoundaryCondition> FlowBoundaryConditions;
    typedef std::vector<SaturationBoundaryCondition> SaturationBoundaryConditions;


} // namespace Dune


#endif // OPENRS_BOUNDARYCONDITIONS_HEADER
