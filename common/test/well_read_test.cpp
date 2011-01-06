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
    if (parser.hasField("COMPDAT")) {
	parser.getCOMPDAT().write(std::cout);    
    }
    if (parser.hasField("WCONINJE")) {
	parser.getWCONINJE().write(std::cout);    
    }
    if (parser.hasField("WCONPROD")) {
	parser.getWCONPROD().write(std::cout);
    }
    if (parser.hasField("WELTARG")) {
	parser.getWELTARG().write(std::cout);    
    }
    return 0;
}
