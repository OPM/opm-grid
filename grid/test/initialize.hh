#include<dune/grid/common/referenceelements.hh>

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
    GridView gridView = grid.leafView();

    // iterate through leaf grid an evaluate c0 at cell center
    LeafIterator endit = gridView.template end<0>();
    for (LeafIterator it = gridView.template begin<0>(); it!=endit; ++it) {
        /*
        // get geometry type
        Dune::GeometryType gt = it->type();

        // get cell center in reference element
        const Dune::FieldVector<ct,dim>& 
        local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);

        // get global coordinate of cell center
        Dune::FieldVector<ct,dimworld> global = 
        it->geometry().global(local);
        */
	Dune::FieldVector<ct, dimworld> centroid = it->geometry().position();

        // initialize cell concentration
        c[mapper.map(*it)] = c0(centroid);
    }
}
