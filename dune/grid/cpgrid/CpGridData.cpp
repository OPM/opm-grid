#include"config.h"
#include"CpGridData.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"
#include"OrientedEntityTable.hpp"
#include"Indexsets.hpp"
#include"PartitionTypeIndicator.hpp"
#include <dune/grid/common/GridPartitioning.hpp>
#include <dune/common/parallel/remoteindices.hh>
#include <algorithm>

namespace Dune
{
namespace cpgrid
{


CpGridData::CpGridData(const CpGridData& g)
  : partition_type_indicator_(new PartitionTypeIndicator(*this)), ccobj_(g.ccobj_)
{}

CpGridData::CpGridData()
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{}

#if HAVE_MPI
CpGridData::CpGridData(MPI_Comm comm)
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(comm), use_unique_boundary_ids_(false)
{}
#endif

CpGridData::CpGridData(CpGrid& grid)
  : index_set_(new IndexSet(*this)),   local_id_set_(new IdSet(*this)),
    global_id_set_(new GlobalIdSet(local_id_set_)),  partition_type_indicator_(new PartitionTypeIndicator(*this)),
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

template<class E>
struct AttributeDataHandle
{
    typedef std::pair<std::size_t,char> DataType;
    
    template<class T>
    AttributeDataHandle(int rank, const std::vector<char>& indicator,
                        std::vector<std::map<std::size_t, char> >& vals,
                        const T& cell_to_entity)
        : rank_(rank), indicator_(indicator), vals_(vals),
          c2e_(cell_to_entity)
    {}
    bool fixedSize()
    {
        return true;
    }
    std::size_t size(std::size_t i)
    {
        return c2f_[i].size();
    }
    template<class B>
    void gather(B buffer, std::size_t i)
    {
        typedef typename SparseTable<E>::row_type::const_iterator RowIter;
        for(RowIter f=c2f[i].begin(), fend=c2f[i].end();
            f!=fend; ++f)
        {
            buffer.write(std::make_pair(rank_,indicator_[f->index()]));
        }
    }

    template<class B>
    void scatter(B buffer, std::size_t i, std::size_t s)
    {
        typedef typename SparseTable<E>::row_type::const_iterator RowIter;
        for(RowIter f=c2f[i].begin(), fend=c2f[i].end();
            f!=fend; ++f, --s)
        {
            std::pair<int,char> rank_attr;
            buffer.read(rank_attr);
            vals[f->index()].insert(rank_attr);
        }
    }
    int rank_;
    const std::vector<char>& indicator;
    std::vector<std::map<std::size_t, char> >& vals;
    const SparseTable<E>& c2f_;
};
 

template<class T, class Functor, class FromSet, class ToSet>
struct InterfaceFunctor
{
    InterfaceFunctor(std::map<int,std::pair<T,T> >& m)
        : map_(m)
    {}
    void operator()(int rank, std::size_t index PartitionType mine, PartitionType other)
    {
        if(from.contains(mine) && to.contains(other))
            func(map_[rank]->first, index);
        if(from.contains(other) && to.contains(mine))
            func(map_[rank]->second, index);
    }
    std::map<int,std::pair<T,T>& map_;
    FromSet from;
    ToSet   to;
    Functor func;
};

struct InterfaceIncrementor
{
    template<class T>
    void operator()(T& t, std::size_t)
    {
        ++t;
    }    
};

struct InterfaceAdder
{
    template<class T>
    void operator()(InterfaceInformation& inf, std::size_t index)
    {
        info->add(index);
    }
};

template<class T0, class T1, class T2, class T3, class T4>
struct InterfaceTupleFunctor<tuple<T0, T1, T2, T3, T4> >
{
    typedef tuple< T0, T1, T2, T3, T4 > MyTuple;
    InterfaceTupleFunctor(MyTuple& t)
    : t_(&t)
    {}
    
    void operator()(int rank, std::size_t index, PartitionType mine, PartitionType other)
    {
        get<0>(t)(rank, index, mine, other);
        get<1>(t)(rank, index, mine, other);
        get<2>(t)(rank, index, mine, other);
        get<3>(t)(rank, index, mine, other);
        get<4>(t)(rank, index, mine, other);
    }
};

struct Converter
{
    typedef EnumItem<PartitionType, InteriorEntity> Interior;
    typedef EnumItem<PartitionType, BorderEntity> Border;
    typedef EnumItem<PartitionType, OverlapEntity> Overlap;
    typedef EnumItem<PartitionType, Front> Front;
    typedef Combine<Interior,Border> InteriorBorder;
    typedef Combine<Overlap,Front> OverlapFront;

    typedef SourceTuple     <InteriorBorder, InteriorBorder, Overlap     , Overlap, AllSet>
    typedef DestinationTuple<InteriorBorder, AllSet,         OverlapFront, AllSet,  AllSet>
};

/**
 * \brief Reserves space for an interface.
 * \tparam i The index of the interface.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
template<std::size_t i>
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
             interfaces)
{
    get<i>(interfaces)->first.reserve(size[i].first);
    get<i>(interfaces)->second.reserve(size[i].second);
}

/**
 * \brief Reserves space for the interfaces.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
             interfaces)
{
    reserve<0>(sizes, interfaces);
    reserve<1>(sizes, interfaces);
    reserve<2>(sizes, interfaces);
    reserve<3>(sizes, interfaces);
    reserve<4>(sizes, interfaces);
}

/**
 * \brief A functor that calculates the size of the interface.
 * \tparam i The indentifier of the interface.
 */
template<std::size_t i>
struct SizeFunctor : 
    public InterfaceFunctor<std::size_t, InterfaceIncrementor,
                            tuple_element<i,Converter::SourceTuple>::type,
                            tuple_element<i,Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<std::size_t, InterfaceIncrementor,
                             tuple_element<i,Converter::SourceTuple>::type,
                             tuple_element<i,Converter::DestinationTuple>::type>
    Base;
    SizeFunctor(std::map<int,std::pair<std::size_t,std::size_t> >& m)
        :Base(m)
    {}
};

/**
 * \brief A functor that adds indices to the interface.
 * \tparam i The indentifier of the interface.
 */
template<std::size_t i>
struct AddFunctor :
    public InterfaceFunctor<InterfaceInformation, InterfaceAdder,
                            tuple_element<i,Converter::SourceTuple>::type,
                            tuple_element<i,Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<InterfaceInformation, InterfaceAdder,
                             tuple_element<i,Converter::SourceTuple>::type,
                             tuple_element<i,Converter::DestinationTuple>::type>
    Base;
    SizeFunctor(std::map<int,std::pair<InterfaceInformation,InterfaceInformation> >& m)
        : Base(m)
    {}
};

/**
 * \brief Applies a functor the each pair of the interface.
 * \tparam Functor The type of the functor to apply.
 * \param attributes[in] A vector that contains for each index a map from other
 * process ranks to the attribute there.
 * \param my_attributes[in] A vector with the attributes of each index on this process.
 * \param func The functor.
 */
template<class Functor>
void iterate_over_attributes(std::vector<std::map<int,char> >& attributes,
                             std::vector<char>& my_attributes, Functor& func)
{
    typedef typename std::vector<std::map<int,char> >::const_iterator Iter;
    typename std::vector<char>::const_iterator mine=my_attributes.begin();
    for(Iter begin=attributes.begin(), i=begin, end=attributes.end(); i!=end; ++i, ++mine)
    {
        typedef typename std::map<int,char>::const_iterator MIter;
        for(MIter m=i->begin(), mend=m.end(); m!=mend; ++m)
            func(m->first,i-begin,*mine, m->second);
    }
}

        
/**
 * \brief Creates the communication interface for either faces or points.
 * \param attributes[in] A vector that contains for each index a map from other
 * process ranks to the attribute there.
 * \param my_attributes[in] A vector with the attributes of each index on this process.
 * \param[out] interfaces The tuple with the interface maps for communication.
 */
void createInterfaces(std::vector<std::map<int,char> >& attributes,
                      std::vector<char>& my_attributes,
                      tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                      interfaces)
{
    // calculate sizes
    std::vector<std::map<int,std::pair<std::size_t,std::size_t> > > sizes(5);
    typedef tuple<SizeFunctor<0>,SizeFunctor<1>,SizeFunctor<2>,SizeFunctor<3>,
                  SizeFunctor<4> > SizeTuple;
    SizeTuple
        size_functor_tuple(SizeFunctor<0>(sizes[0]),
                           SizeFunctor<1>(sizes[1]),
                           SizeFunctor<2>(sizes[2]),
                           SizeFunctor<3>(sizes[3]),
                           SizeFunctor<4>(sizes[4]));
    InterfaceTupleFunctor<SizeTuple> size_functor(size_functor_tuple);
    iterate_over_attributes(attributes, my_attributes, size_functor);
    // reserve space
    reserve(sizes, interfaces);
    // add indices to the interface
    typedef tuple<AddFunctor<0>,AddFunctor<1>,AddFunctor<2>,AddFunctor<3>,
                  AddFunctor<4> > AddTuple;
    AddTuple
        add_functor_tuple(AddFunctor<0>(get<0>(interfaces)),
                          AddFunctor<1>(get<1>(interfaces)),
                          AddFunctor<2>(get<2>(interfaces)),
                          AddFunctor<3>(get<3>(interfaces)),
                          AddFunctor<4>(get<4>(interfaces)));
    InterfaceTupleFunctor<AddTuple> add_functor(add_functor_tuple);
    iterate_over_attributes(attributes, my_attributes, add_functor);
                                
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
        ParallelIndexSet* indexset;
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
    cell_indexset = new ParallelIndexSet;
    cell_counter.indexset=indexset;
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
        // Therfore we use an ugly cast to base class here.
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
                // The statement below results in all faces having two neighbours
                // except for those at the domain boundary.
                // Note that along the front partition there are invalid neighbours
                // marked with index std::numeric_limits<int>::max()
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
    // Compute the partition type for cell
    partition_type_indicator_->cell_indicator_.resize(indexset.size());
    for(IndexSet::const_iterator i=indexset.begin(), end=indexset.end();
            i!=end; ++i)
    {
        partition_type_indicator_->cell_indicator_[i->local()]=
            i->local().attribute()==OwnerOverlapCopyAttributeSet::owner?
            InteriorEntity:OverlapEntity;
    }
    
    // Compute partition type for points
    // We initialize all point with interior. Then we loop over the faces. If a face is of 
    // type border, then the type of the point is overwritten with border. In the other cases
    // we set type of the point to the one of the face as long as the type of the point is 
    // not boder.
    partition_type_indicator_->point_indicator_.resize(geometry_.geomVector<3>().size(),
                                                      InteriorEntity);
    for(int i=0; i<face_to_point_.size(); ++i)
    {
        for(auto p=view_data.face_to_point_[i].begin(), 
                pend=view_data.face_to_point_[i].end(); p!=pend; ++p)
        {
            if(partition_type_indicator_->point_indicator_[*p]!=BorderEntity)
                partition_type_indicator_->point_indicator_[*p]=
                    partition_type_indicator_->getFacePartitionType(i);
        }
    }

    // Compute the interface information for cells
    Dune::RemoteIndices<IndexSet> cell_remote_indices.setIndexSets(cell_indexset, cellindexset, ccobj_);
    get<InteriorBorder_All_Interface>(cell_interfaces)
        .build(cell_remote_indices, EnumItem<AttributeSet, AttributeSet::owner>(),
               AllSet<AttributeSet>());
    get<Overlap_OverlapFront_Interface>(cell_interfaces)
        .build(cell_remote_indices, EnumItem<AttributeSet, AttributeSet::overlap>(),
               EnumItem<AttributeSet, AttributeSet::overlap>());
     get<Overlap_All_Interface>(cell_interfaces)
        .build(cell_remote_indices, EnumItem<AttributeSet, AttributeSet::overlap>(),
                                 AllSet<AttributeSet>());
    get<All_All_Interface>(cell_interfaces)
        .build(cell_remote_indices, AllSet<AttributeSet>(), AllSet<AttributeSet>());

    // Now we use the all_all communication of the cells to compute which faces and points
    // are also present on other processes and with what attribute.
    std::vector<std::map<int,char> > face_attributes(noExistingFaces);
    FaceAttributeDataHandle face_handle(ccobj_.rank(), 
                                        partition_type_indicator_->face_indicator,
                                        face_attributes,
                                        cell_to_face_);
    Dune::VariableSizeCommunicator<> comm(get<All_All_Interface>(cell_interfaces));
    comm.forward(face_handle);
    createInterfaces(face_attributes, face_indicator, face_interfaces);
    std::vector<std::map<int,char> >().swap(face_attributes);
    std::vector<std::map<int,char> > point_attributes(noExistingPoints);
    PointAttributeDataHandle point_handle(ccobj_.rank(), 
                                          partition_type_indicator_->point_indicator_,
                                          point_attributes,
                                          cell_to_point_);
    comm.forward(point_handle);
    createInterfaces(point_attributes, partition_type_indicator_->point_indicator_,
                     point_interfaces);
    
    
}

} // end namespace cpgrid
} // end namespace Dune
