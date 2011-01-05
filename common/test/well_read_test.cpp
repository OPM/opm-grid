//===========================================================================
//                                                                           
// File: well_read_test.cpp                                                  
//                                                                           
//===========================================================================

#include "../EclipseGridParser.hpp"
#include <iostream>

int main()
{
    const std::string filename = "SPE9.DATA";
    Dune::EclipseGridParser parser(filename);
    if (parser.hasField("WELSPECS")) {
	parser.getWELSPECS().write(std::cout);    
    }
    return 0;
}
