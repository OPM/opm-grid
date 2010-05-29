// $Id: test_diffusion.cc 2173 2009-06-19 11:58:37Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Markus Wolff                 *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/CpGrid.hpp>

template<int dim>
struct ElementLayout
{
    bool contains (Dune::GeometryType gt)
    {
	// std::cout << "dim = " << dim << "    gt.dim() = " << gt.dim() << "     gt = " << gt << std::endl;
        return dim == gt.dim();//gt.isSingular() && gt.dim() == 3;
    }
};


int main(int /*argc*/, char** /*argv*/)
{
    try {
        typedef Dune::CpGrid Grid;

        int refinement = 1;
        Grid grid;
        Dune::array<int   , 3> dims;
        std::fill(dims.begin(), dims.end(), 1 << refinement);
        Dune::array<double, 3> cell_sz;
        std::fill(cell_sz.begin(), cell_sz.end(), 1.0 / (1 << refinement));
        grid.createCartesian(dims, cell_sz);
        
        typedef Grid::LeafGridView GridView;
        GridView gridView(grid.leafView());

        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> ElementMapper;
        ElementMapper elementMapper(gridView);

        typedef GridView::IndexSet IndexSet;
        const IndexSet& indexSet(gridView.indexSet());
        
        typedef GridView::Codim<0>::Iterator ElementIterator;
        ElementIterator endEIt = gridView.end<0>();
        for (ElementIterator eIt = gridView.begin<0>(); eIt != endEIt; ++eIt) {
            int mapperIdx = elementMapper.map(*eIt);
            int indexSetIdx = indexSet.index(*eIt);
	    if (mapperIdx != indexSetIdx) {
		std::cerr << "Mismatched mapper and indexset indices: " << mapperIdx << " vs. " << indexSetIdx << '\n';
		return EXIT_FAILURE;
	    }
	    // std::cout << "mapperIdx = " << mapperIdx << ", indexSetIdx = " << indexSetIdx << std::endl;
        }
	return EXIT_SUCCESS;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
    return EXIT_FAILURE;
}

