/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

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
#ifndef OPM_GRID_ERROR_MACROS_HPP
#define OPM_GRID_ERROR_MACROS_HPP

#if HAVE_OPM_COMMON
#include <opm/common/ErrorMacros.hpp>
#else

#include <iostream>
#include <string>

#define OPM_THROW(Exception, message)                    \
    do {                                                       \
        std::string oss_ = std::string{"["} + __FILE__ + ":" + \
                           std::to_string(__LINE__) + "] " +   \
                           message;                            \
        throw Exception(oss_);                                 \
    } while (false)

#define OPM_THROW_NOLOG OPM_THROW

// throw an exception if a condition is true
#define OPM_ERROR_IF(condition, message) do {if(condition){ OPM_THROW(std::logic_error, message);}} while(false)

#define OPM_MESSAGE(x) do { std::cerr << x << std::endl; } while(false)

#endif

#endif // OPM_GRID_ERROR_MACROS_HPP
