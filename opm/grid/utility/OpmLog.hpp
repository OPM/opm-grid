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
#ifndef OPM_GRID_LOG_HPP
#define OPM_GRID_LOG_HPP

#if HAVE_OPM_COMMON
#include <opm/common/OpmLog/OpmLog.hpp>
#else

#include <iostream>
#include <string_view>

namespace Opm::OpmLog {
  static inline void error(std::string_view msg)
  {
      std::cerr << msg << std::endl;
  }

  static inline void info(std::string_view msg)
  {
      std::cout << msg << std::endl;
  }

  static inline void warning(std::string_view msg)
  {
      std::cout << msg << std::endl;
  }

  static inline void warning(std::string_view tag, std::string_view msg)
  {
      std::cout << '[' << tag << "]: " << msg << std::endl;
  }

  static inline void setupSimpleDefaultLogging() {}
}
#endif

#endif // OPM_GRID_LOG_HPP
