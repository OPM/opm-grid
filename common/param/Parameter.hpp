//===========================================================================
//
// File: Parameter.hpp
//
// Created: Tue Jun  2 16:00:21 2009
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

#ifndef OPENRS_PARAMETER_HEADER
#define OPENRS_PARAMETER_HEADER

#include <string>
#include <sstream>

#include <dune/common/param/ParameterMapItem.hpp>
#include <dune/common/param/ParameterStrings.hpp>

namespace Dune {
    /// See ParameterGroup.hpp for how to use the parameter system
    namespace parameter {
	class Parameter : public ParameterMapItem {
	public:
	    virtual ~Parameter() {}
	    virtual std::string getTag() const {return ID_xmltag__param;}
	    Parameter(const std::string& value, const std::string& type)
                : value_(value), type_(type) {}
	    std::string getValue() const {return value_;}
	    std::string getType() const {return type_;}
	private:
	    std::string value_;
	    std::string type_;
	};

	std::string correct_parameter_tag(const ParameterMapItem& item);
	std::string correct_type(const Parameter& parameter,
                                 const std::string& type);

	template<>
	struct ParameterMapItemTrait<int> {
	    static int convert(const ParameterMapItem& item,
                               std::string& conversion_error)
            {
		conversion_error = correct_parameter_tag(item);
		if (conversion_error != "") {
		    return 0;
		}
		const Parameter& parameter = dynamic_cast<const Parameter&>(item);
		conversion_error = correct_type(parameter, ID_param_type__int);
		if (conversion_error != "") {
		    return 0;
		}
		std::stringstream stream;
		stream << parameter.getValue();
		int value;
		stream >> value;
		if (stream.fail()) {
		    conversion_error = "Conversion to '" +
                                       ID_param_type__int +
                                       "' failed. Data was '" +
                                       parameter.getValue() + "'.\n";
		    return 0;
		}
		return value;
	    }
	    static std::string type() {return ID_param_type__int;}
	};

	template<>
	struct ParameterMapItemTrait<double> {
	    static double convert(const ParameterMapItem& item,
                                  std::string& conversion_error)
            {
		conversion_error = correct_parameter_tag(item);
		if (conversion_error != "") {
		    return 0.0;
		}
		const Parameter& parameter = dynamic_cast<const Parameter&>(item);
		conversion_error = correct_type(parameter, ID_param_type__float);
		if (conversion_error != "") {
		    return 0.0;
		}
		std::stringstream stream;
		stream << parameter.getValue();
		double value;
		stream >> value;
		if (stream.fail()) {
		    conversion_error = "Conversion to '" +
                                       ID_param_type__float +
                                       "' failed. Data was '" +
                                       parameter.getValue() + "'.\n";
		    return 0.0;
		}
		return value;
	    }
	    static std::string type() {return ID_param_type__float;}
	};

	template<>
	struct ParameterMapItemTrait<bool> {
	    static bool convert(const ParameterMapItem& item,
                                std::string& conversion_error)
            {
		conversion_error = correct_parameter_tag(item);
		if (conversion_error != "") {
		    return false;
		}
		const Parameter& parameter = dynamic_cast<const Parameter&>(item);
		conversion_error = correct_type(parameter, ID_param_type__bool);
		if (conversion_error != "") {
		    return false;
		}
		if (parameter.getValue() == ID_true) {
		    return true;
		} else if (parameter.getValue() == ID_false) {
		    return false;
		} else {
		    conversion_error = "Conversion failed. Data was '" +
                                       parameter.getValue() +
                                       "', but should be one of '" +
                                       ID_true + "' or '" + ID_false + "'.\n";
		    return false;
		}
	    }
	    static std::string type() {return ID_param_type__bool;}
	};

	template<>
	struct ParameterMapItemTrait<std::string> {
	    static std::string convert(const ParameterMapItem& item,
                                       std::string& conversion_error)
            {
		conversion_error = correct_parameter_tag(item);
		if (conversion_error != "") {
		    return "";
		}
		const Parameter& parameter = dynamic_cast<const Parameter&>(item);
		conversion_error = correct_type(parameter, ID_param_type__string);
		if (conversion_error != "") {
		    return "";
		}
		return parameter.getValue();
	    }
	    static std::string type() {return ID_param_type__string;}
	};
    } // namespace parameter
} // namespace Dune
#endif  // OPENRS_PARAMETER_HPP
