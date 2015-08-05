#include <dune/common/version.hh>
#include<dune/geometry/referenceelements.hh>
#include "transportproblem2.hh"

//! initialize the vector of unknowns with initial value
template<class G, class M, class V>
void initialize(const G& grid, const M& mapper, V& c)
{
    // first we extract the dimensions of the grid
    const int dim = G::dimension;
    const int dimworld = G::dimensionworld;

    // type used for coordinates in the grid
    typedef typename G::ctype ct;

    // type of grid view on leaf part
    typedef typename G::LeafGridView GridView;

    // leaf iterator type
    typedef typename GridView::template Codim<0>::Iterator LeafIterator;

    // get grid view on leaf part
    GridView gridView = grid.leafGridView();

    // iterate through leaf grid an evaluate c0 at cell center
    LeafIterator endit = gridView.template end<0>();
    for (LeafIterator it = gridView.template begin<0>(); it!=endit; ++it) {
        // get geometry type
        Dune::GeometryType gt = it->type();

        // get cell center in reference element
        const Dune::FieldVector<ct,dim>&
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
        local = Dune :: ReferenceElements<ct,dim >::general(gt).position(0,0);
#else
        local = Dune :: GenericReferenceElements<ct,dim >::general(gt).position(0,0);
#endif

        // get global coordinate of cell center
        Dune::FieldVector<ct,dimworld> global =
            it->geometry().global(local);
        // Dune::FieldVector<ct, dimworld> global = it->geometry().position();

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
        // initialize cell concentration
        c[mapper.index(*it)] = c0(global);
#else
        // initialize cell concentration
        c[mapper.map(*it)] = c0(global);
#endif
    }
}
