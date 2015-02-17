// Warning suppression for Dune includes.
#include <opm/core/utility/platform_dependent/disable_warnings.h>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/core/utility/platform_dependent/reenable_warnings.h>

#include <stdio.h>

template<class G, class V>
void vtkout (const G& grid, const V& c, const char* name, int k, double time=0.0, int rank=0)
{
  Dune::VTKWriter<typename G::LeafGridView> vtkwriter(grid.leafGridView());
  char fname[128];
  char sername[128];
  sprintf(fname,"%s-%05d",name,k);
  sprintf(sername,"%s.series",name);
  vtkwriter.addCellData(c,"celldata");
  vtkwriter.write(fname,Dune::VTK::ascii);
  if ( rank == 0) 
  {
    std::ofstream serstream(sername, (k==0 ? std::ios_base::out : std::ios_base::app));
    serstream << k << " " << fname << ".vtu " << time << std::endl;
    serstream.close();
  }
}
