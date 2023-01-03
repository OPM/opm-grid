#include "DataHandleWrappers.hpp"

#include <iostream>

bool Dune::cpgrid::PointViaCellWarner::printWarn = true;

void Dune::cpgrid::PointViaCellWarner::warn()
{
    if (printWarn)
    {
        std::cerr << "Communication of variable data attached to points is "
                  << "not fully supported. Your code/handle must not use the "
                  << "last parameter of "
                  << "DataHandle::scatter(B& buffer, E& entity, int size) "
                  << "as it will not be correct!\n";
        printWarn = false;
    }
}
