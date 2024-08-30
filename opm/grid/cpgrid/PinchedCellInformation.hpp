/*
  Copyright 2024 Statoil ASA

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
#ifndef OPM_PINCHEDCELLINFORMATION_HEADER_INCLUDED
#define OPM_PINCHEDCELLINFORMATION_HEADER_INCLUDED
#include<array>
#include<tuple>

namespace Opm
{

class PinchedCellInformation
{
public:
    using CenterContainer = std::vector<std::pair<std::array<double,3>,std::array<double,3>>>;
    using iterator = CenterContainer::const_iterator;
    void push_back(const std::pair<std::array<double,3>,std::array<double,3>>& centers)
    {
        cellAndBottomFaceCenter_.push_back(centers);
    }
    iterator begin() const
    {
        return cellAndBottomFaceCenter_.begin();
    }
    iterator end() const
    {
        return cellAndBottomFaceCenter_.end();
    }
    
private:
    CenterContainer cellAndBottomFaceCenter_;
};
}
#endif
