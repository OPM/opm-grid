#include"config.h"
#include"CpGridData.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"
#include"OrientedEntityTable.hpp"
#include"Indexsets.hpp"
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/grid/common/GridPartitioning.hpp>
#include <algorithm>

namespace Dune
{
namespace cpgrid
{
CpGridData::CpGridData()
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)),
      ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{}

#if HAVE_MPI
CpGridData::CpGridData(MPI_Comm comm)
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)),
      ccobj_(comm), use_unique_boundary_ids_(false)
{}
#endif

CpGridData::CpGridData(CpGrid& grid)
  : index_set_(new IndexSet(*this)),   local_id_set_(new IdSet(*this)),
    global_id_set_(new GlobalIdSet(local_id_set_)),
    ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{}
CpGridData::~CpGridData()
{
    delete index_set_;
    delete local_id_set_;
    delete global_id_set_;
}

void CpGridData::computeUniqueBoundaryIds()
{
    // Perhaps we should make available a more comprehensive interface
    // for EntityVariable, so that we don't have to build a separate
    // vector and assign() to unique_boundary_ids_ at the end.
    int num_faces = face_to_cell_.size();
    std::vector<int> ids(num_faces, 0);
    int count = 0;
    for (int i = 0; i < num_faces; ++i) {
        cpgrid::EntityRep<1> face(i, true);
        if (face_to_cell_[face].size() == 1) {
            // It's on the boundary.
            // Important! Since boundary ids run from 1 to n,
            // we use preincrement instead of postincrement below.
            ids[i] = ++count;
        }
    }
    unique_boundary_ids_.assign(ids.begin(), ids.end());
#ifdef VERBOSE
    std::cout << "computeUniqueBoundaryIds() gave all boundary intersections\n"
              << "unique boundaryId()s ranging from 1 to " << count << std::endl;
#endif
}

int CpGridData::size(int codim) const
{
    switch (codim) {
    case 0: return cell_to_face_.size();
    case 1: return 0;
    case 2: return 0;
    case 3: return geomVector<3>().size();
    default: return 0;
    }
}


 // A functor that counts existent entries and renumbers them.
struct CountExistent
{
    CountExistent() : count() {}

    void operator()(int& i)
    {
        if(i < std::numeric_limits<int>::max())
            count++;
    }
    int count;
};

struct AssignAndIncrement
{
    AssignAndIncrement() : i_(){}
    void operator()(int& val){ if(val<std::numeric_limits<int>::max()) val=i_++; }
    int i_;
} assigner;

/** 
 * @brief Counts the number of global ids and sets them up.
 * @param indicator A vector indicating whether an entity exists.
 * @param ids A vector to the the global ids in.
 * @param idSet The idSet of the global grid.
 * @return the number of entities that exist.
 */
template<int codim>
int setupAndCountGlobalIds(const std::vector<int>& indicator, std::vector<int>& ids,
                           const IdSet& idSet)
{
    int count = std::count_if(indicator.begin(),
                              indicator.end(), 
                              std::bind2nd(std::less<int>(),
                                           std::numeric_limits<int>::max()));
    ids.resize(count);
    typedef typename std::vector<int>::const_iterator VIter;
    for(VIter ibegin=indicator.begin(), i=ibegin, iend= indicator.end(); 
        i!=iend; ++i)
    {
        ids[*i]=idSet.id(EntityRep<codim>(i-ibegin,true));
    }
    return count;
}

void CpGridData::distributeGlobalGrid(const CpGrid& grid,
                                      const CpGridData& view_data,
                                      const std::vector<int>& cell_part)
{
    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>& ccobj=ccobj_;
    std::unique_ptr<CpGridData> grid_data(new CpGridData);
    int my_rank=ccobj.rank();
#ifdef DEBUG
    int size=ccobj.size();
    typedef std::vector<int>::const_iterator Iterator;
    for(Iterator i=cell_part.begin(); i!= cell_part.end(); ++i)
        if(*i>=size)
            OPM_THROW(std::runtime_error, "rank for cell is too big");
#endif
    // vector with the set of ranks that
    std::vector<std::set<int> > overlap; 
    
    //Dune::OwnerOverlapCopyCommunication<int,Dune::ParallelLocalIndex> ownerOverlap;
    typedef OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
    typedef Dune::ParallelIndexSet<int,Dune::ParallelLocalIndex<AttributeSet> > IndexSet;
    overlap.resize(cell_part.size());
    addOverlapLayer(grid, cell_part, overlap, my_rank, false);
    // count number of cells
    struct CellCounter
    {
        int myrank;
        int count;
        typedef std::set<int>::const_iterator NeighborIterator;
        std::set<int> neighbors;
        typedef Dune::ParallelLocalIndex<AttributeSet> Index;
        IndexSet* indexset;
        std::vector<int> global2local;
        /**
         * @brief Adds an index with flag owner to the index set.
         *
         * An index is added if the rank is our rank.
         * @param i The global index
         * @param rank The rank that owns the index.
         */
        void operator() (int i, int rank)
        {
            if(rank==myrank)
            {
                global2local.push_back(count);
                indexset->add(i, Index(count++, OwnerOverlapCopyAttributeSet::owner, true));
            }
            else
                global2local.push_back(std::numeric_limits<int>::max());
        }
        /**
         * @brief Adds an index with flag overlap to the index set.
         * @param i The global index.
         * @param ov The set of ranks where the index is in the overlap
         * region.
         * @param rank The rank that owns this index.
         */
        void operator() (int i, const std::set<int>& ov, int rank)
        {
            std::set<int>::const_iterator iter=
                ov.find(myrank);
            if(iter!=ov.end())
            {
                global2local.push_back(count);
                indexset->add(i, Index(count++, OwnerOverlapCopyAttributeSet::overlap, true));
                neighbors.insert(rank);
            }
            else
                global2local.push_back(std::numeric_limits<int>::max());
            
        }
    } cell_counter;
    cell_counter.myrank=my_rank;
    cell_counter.count=0;
    cell_counter.global2local.reserve(cell_part.size());
    IndexSet indexset;
    cell_counter.indexset=&indexset;
    // set up the index set.
    cell_counter.indexset->beginResize();
    typedef std::vector<std::set<int> >::const_iterator OIterator;
    std::vector<int>::const_iterator ci=cell_part.begin();
    for(OIterator end=overlap.end(), begin=overlap.begin(), i=begin; i!=end; ++i, ++ci)
    {
        if(i->size())
            cell_counter(i-begin, *i, *ci);
        else
            cell_counter(i-begin, *ci);
    }
    cell_counter.indexset->endResize();
    // setup the remote indices.
    //typedef Dune::OwnerOverlapCopyCommunication<int,Dune::ParallelLocalIndex>::RemoteIndices RemoteIndices;
    typedef Dune::RemoteIndices<IndexSet> RemoteIndices;
    typedef RemoteIndexListModifier<RemoteIndices::ParallelIndexSet, RemoteIndices::Allocator,
                                    false> Modifier;
    typedef RemoteIndices::RemoteIndex RemoteIndex;
    RemoteIndices remoteIndices;

    // Create a map of ListModifiers
    { //extra scope to call destructor of the Modifiers
        std::map<int,Modifier> modifiers;
        for(CellCounter::NeighborIterator n=cell_counter.neighbors.begin(), end=cell_counter.neighbors.end();
            n != end; ++n)
            modifiers.insert(std::make_pair(*n, remoteIndices.getModifier<false,false>(*n)));
        // Insert remote indices. For each entry in the index set, see wether there are overlap occurences and add them.
        for(IndexSet::const_iterator i=indexset.begin(), end=indexset.end();
            i!=end; ++i)
        {
            std::set<int>::const_iterator iter=overlap[i->local()].find(my_rank);
            if(iter!=overlap[i->local()].end())
                modifiers[cell_part[i->local()]]
                    .insert(RemoteIndex(OwnerOverlapCopyAttributeSet::owner,&(*i)));
        }
    }

    
    // We can identify existing cells with the help of the index set.
    // Now we need to compute the existing faces and points. Either exist
    // if they are reachable from an existing cell.
    // We use std::numeric_limits<int>::max() to indicate non-existent entities.
    std::vector<int> face_indicator(view_data.geometry_.geomVector<1>().size(),
                                    std::numeric_limits<int>::max());
    std::vector<int> point_indicator(view_data.geometry_.geomVector<3>().size(), 
                                     std::numeric_limits<int>::max());
    for(IndexSet::iterator i=indexset.begin(), end=indexset.end();
            i!=end; ++i)
    {
        typedef boost::iterator_range<const EntityRep<1>*>::iterator RowIter;
        int row_index=i->global();
        // Somehow g++-4.4 does not find functions of father even if we
        // change inheritance of OrientedEntityTable to public.
        // Therfore we use an ugly task to base class here.
        const Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        for(RowIter f=c2f[row_index].begin(), fend=c2f[row_index].end();
            f!=fend; ++f)
        {
            int findex=f->index();
            --face_indicator[findex];
            // points reachable from a cell exist, too.
            for(auto p=view_data.face_to_point_[findex].begin(), 
                    pend=view_data.face_to_point_[findex].end(); p!=pend; ++p)
            {
                assert(p>=0);
                --point_indicator[*p];
            }
        } 
    }
    
       
    // renumber  face and point indicators
    std::for_each(face_indicator.begin(), face_indicator.end(), AssignAndIncrement());
    std::for_each(point_indicator.begin(), point_indicator.end(), AssignAndIncrement());

    std::vector<int> map2GlobalFaceId;
    int noExistingFaces = setupAndCountGlobalIds<1>(face_indicator, map2GlobalFaceId,
                                                    *view_data.local_id_set_);
    std::vector<int> map2GlobalPointId;
    int noExistingPoints = setupAndCountGlobalIds<3>(point_indicator, map2GlobalPointId,
                                                     *view_data.local_id_set_);
    std::vector<int> map2GlobalCellId(indexset.size());
    for(IndexSet::const_iterator i=indexset.begin(), end=indexset.end();
        i!=end; ++i)
    {
        map2GlobalCellId[i->local()]=view_data.local_id_set_->id(EntityRep<0>(i->global(), true));
    }

    global_id_set_->swap(map2GlobalCellId, map2GlobalFaceId, map2GlobalPointId);
    // Set up the new topology arrays
    EntityVariable<cpgrid::Geometry<3, 3>, 0> cell_geom;
    std::vector<cpgrid::Geometry<3, 3> > tmp_cell_geom(indexset.size());
    auto global_cell_geom=view_data.geomVector<0>();
    std::vector<int> global_cell;
    global_cell.resize(indexset.size());

    // Copy the existing cells.
    for(auto i=indexset.begin(), end=indexset.end(); i!=end; ++i)
    {
        tmp_cell_geom[i->local()]=static_cast<std::vector<cpgrid::Geometry<3, 3> >&>(global_cell_geom)[i->global()];
        global_cell[i->local()]=view_data.global_cell_[i->global()];
    }
    static_cast<std::vector<cpgrid::Geometry<3, 3> >&>(cell_geom).swap(tmp_cell_geom);

    // count the existing faces, renumber, and allocate space.
    EntityVariable<cpgrid::Geometry<2, 3>, 1> face_geom;
    std::vector<cpgrid::Geometry<2, 3> > tmp_face_geom(noExistingFaces);
    std::vector<enum face_tag> tmp_face_tag(noExistingFaces);
    std::vector<PointType> tmp_face_normals(noExistingFaces);
    const std::vector<cpgrid::Geometry<2, 3> >& global_face_geom=view_data.geomVector<1>();
    const std::vector<enum face_tag>& global_face_tag=view_data.face_tag_;
    const std::vector<PointType>& global_face_normals=view_data.face_normals_;
    
    // Now copy the face geometries that do exist.
    auto f = tmp_face_geom.begin();
    auto ft = tmp_face_tag.begin();
    auto fn = tmp_face_normals.begin();
    for(auto begin=face_indicator.begin(), fi=begin,  end=face_indicator.end(); fi!=end; 
        ++fi)
    {
        if(*fi<std::numeric_limits<int>::max())
        {
            *f=global_face_geom[fi-begin];
            *ft=global_face_tag[fi-begin];
            *fn=global_face_normals[fi-begin];
            ++f; ++ft; ++fn;
        }
    }
    static_cast<std::vector<PointType>&>(face_normals_).swap(tmp_face_normals);
    static_cast<std::vector<cpgrid::Geometry<2, 3> >&>(face_geom).swap(tmp_face_geom);
    static_cast<std::vector<enum face_tag>&>(face_tag_).swap(tmp_face_tag);

    // Count the existing points and allocate space
    std::vector<cpgrid::Geometry<0, 3> > tmp_point_geom(noExistingPoints);
    EntityVariable<cpgrid::Geometry<0, 3>, 3> point_geom;
    const std::vector<cpgrid::Geometry<0, 3> >& global_point_geom=view_data.geomVector<3>();
    
    // Now copy the point geometries that do exist.
    auto p= tmp_point_geom.begin();
    for(auto begin=point_indicator.begin(), pi=begin,  end=point_indicator.end(); pi!=end;
        ++pi)
    {
        if(*pi<std::numeric_limits<int>::max())
        {
            *p=global_point_geom[pi-begin];
            ++p;
        }
    }
    // swap the underlying vectors to get data into point_geom
    static_cast<std::vector<cpgrid::Geometry<0, 3> >&>(point_geom).swap(tmp_point_geom);
    
    // Copy the vectors to geometry. There is no other way currently.
    geometry_=cpgrid::DefaultGeometryPolicy(cell_geom, face_geom, 
                                            point_geom);

    // Create the topology information. This is stored in sparse matrix like data structures.
    // First conunt the size of the nonzeros of the cell_to_face data.
    int data_size=0;
    for(auto i=indexset.begin(), end=indexset.end();
        i!=end; ++i)
    {
        const Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        data_size+=c2f.rowSize(i->global());
    }
    
    //- cell_to_face_ : extract owner/overlap rows from cell_to_face_
    // Construct the sparse matrix like data structure.
    //OrientedEntityTable<0, 1> cell_to_face;
    cell_to_face_.reserve(indexset.size(), data_size);
    cell_to_point_.reserve(indexset.size());
    
    for(auto i=indexset.begin(), end=indexset.end();
        i!=end; ++i)
    {
        typedef boost::iterator_range<const EntityRep<1>*>::const_iterator RowIter;
        const Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        auto row=c2f[i->global()];
        // create the new row, i.e. copy orientation and use new face indicator.
        std::vector<EntityRep<1> > new_row(row.size());
        std::vector<EntityRep<1> >::iterator  nface=new_row.begin();
        for(RowIter face=row.begin(), fend=row.end(); face!=fend; ++face, ++nface)
            nface->setValue(face_indicator[face->index()], face->orientation());
        // Append the new row to the matrix
        cell_to_face_.appendRow(new_row.begin(), new_row.end());
        for(int j=0; j<8; ++j)
            cell_to_point_[i->local()][j]=point_indicator[view_data.cell_to_point_[i->global()][j]];
    }

    // Calculate the number of nonzeros needed for the face_to_cell sparse matrix
    // To speed things up, we only use an upper limit here.
    data_size=0;
    const Opm::SparseTable<EntityRep<0> >& f2c=view_data.face_to_cell_;
    for(auto f=face_indicator.begin(), fend=face_indicator.end(); f!=fend; ++f)
        if(*f<std::numeric_limits<int>::max())
            data_size += f2c.rowSize(*f);

    face_to_cell_.reserve(f2c.size(), data_size);
    
    //- face_to cell_ : extract rows that connect to an existent cell
    std::vector<int> cell_indicator(view_data.cell_to_face_.size(),
                                    std::numeric_limits<int>::max());
    for(auto i=indexset.begin(), end=indexset.end(); i!=end; ++i)
        cell_indicator[i->global()]=i->local();
    
    for(auto begin=face_indicator.begin(), f=begin, fend=face_indicator.end(); f!=fend; ++f)
    {
        if(*f<std::numeric_limits<int>::max())
        {
            // face does exist
            std::vector<EntityRep<0> > new_row;
            auto old_row = f2c[f-begin];
            new_row.reserve(old_row.size());
            // push back connected existent cells.
            // for those cells we use the new cell_indicator and copy the 
            // orientation of the old cell.
            for(auto cell = old_row.begin(), cend=old_row.end();cell!=cend; ++cell)
            {
                if(cell_indicator[cell->index()]<std::numeric_limits<int>::max())
                    new_row.push_back(EntityRep<0>(cell_indicator[cell->index()], 
                                                   cell->orientation()));
            }
            face_to_cell_.appendRow(new_row.begin(), new_row.end());
        }
    }
    
    // Compute the number of non zeros of the face_to_point matrix.
    data_size=0;
    for(auto f=face_indicator.begin(), fend=face_indicator.end(); f!=fend; ++f)
        if(*f<std::numeric_limits<int>::max())
            data_size += view_data.face_to_point_.rowSize(*f);

    face_to_point_.reserve(view_data.face_to_point_.size(), data_size);

    //- face_to_point__ : extract row associated with existing faces_
    for(auto begin=face_indicator.begin(), f=begin, fend=face_indicator.end(); f!=fend; ++f)
    {
        if(*f<std::numeric_limits<int>::max())
        {
            // face does exist
            std::vector<int> new_row;
            auto old_row = view_data.face_to_point_[f-begin];
            new_row.reserve(old_row.size());
            for(auto point = old_row.begin(), pend=old_row.end(); point!=pend; ++point)
            {
                assert(point_indicator[*point]<std::numeric_limits<int>::max());
                new_row.push_back(point_indicator[*point]);
            }
            face_to_point_.appendRow(new_row.begin(), new_row.end());
        }
    }
    
    // \TODO - logical_cartesian_size ????
    // - unique_boundary_ids_ : extract the ones that correspond existent faces
    EntityVariable<int, 1> unique_boundary_ids;
    
    if(view_data.unique_boundary_ids_.size())
    {
        // Unique boundary ids are inherited from the global grid.
        auto id=view_data.unique_boundary_ids_.begin();
        unique_boundary_ids.reserve(view_data.face_to_cell_.size());
        for(auto f=face_indicator.begin(), fend=face_indicator.end(); f!=fend; ++f, ++id)
        {
            if(*f<std::numeric_limits<int>::max())
            {
                unique_boundary_ids.push_back(*id);
            }
        }
    }

}

} // end namespace cpgrid
} // end namespace Dune
