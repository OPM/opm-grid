#include"config.h"
#include <algorithm>
#include <map>
#include <vector>
#include"CpGridData.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"
#include"OrientedEntityTable.hpp"
#include"Indexsets.hpp"
#include"PartitionTypeIndicator.hpp"

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <opm/grid/common/GridPartitioning.hpp>
#include <dune/common/parallel/remoteindices.hh>
#include <dune/common/enumset.hh>
#include <opm/grid/utility/SparseTable.hpp>

#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
namespace cpgrid
{


CpGridData::CpGridData(const CpGridData& g)
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)), ccobj_(g.ccobj_)
{
#if HAVE_MPI
    ccobj_=CollectiveCommunication(MPI_COMM_SELF);
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
}

CpGridData::CpGridData()
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{
#if HAVE_MPI
    ccobj_=CollectiveCommunication(MPI_COMM_SELF);
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
}

#if HAVE_MPI
CpGridData::CpGridData(MPI_Comm comm)
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(comm), use_unique_boundary_ids_(false)
{
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
}
#endif

CpGridData::CpGridData(CpGrid&)
  : index_set_(new IndexSet(*this)),   local_id_set_(new IdSet(*this)),
    global_id_set_(new GlobalIdSet(local_id_set_)),  partition_type_indicator_(new PartitionTypeIndicator(*this)),
    ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{
#if HAVE_MPI
    //ccobj_=CollectiveCommunication(MPI_COMM_SELF);
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
}

#if HAVE_MPI
template<class InterfaceMap>
void freeInterfaces(InterfaceMap& map)
{
    typedef typename InterfaceMap::iterator Iter;
    for(Iter i=map.begin(), end=map.end(); i!=end; ++i)
    {
        i->second.first.free();
        i->second.second.free();
    }
}

template<class InterfaceMap>
void freeInterfaces(std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                    interfaces)
{
    freeInterfaces(std::get<0>(interfaces));
    freeInterfaces(std::get<1>(interfaces));
    freeInterfaces(std::get<2>(interfaces));
    freeInterfaces(std::get<3>(interfaces));
    freeInterfaces(std::get<4>(interfaces));
}
#endif

CpGridData::~CpGridData()
{
#if HAVE_MPI
    // code deactivated, because users cannot access face indices and therefore
    // communication on faces makes no sense!
    //freeInterfaces(face_interfaces_);
    freeInterfaces(point_interfaces_);
#endif
    delete index_set_;
    delete local_id_set_;
    delete global_id_set_;
    delete partition_type_indicator_;
}

void CpGridData::populateGlobalCellIndexSet()
{
    cell_indexset_.beginResize();
    for (int index = 0, end = size(0); index != end ; ++index){
        cell_indexset_.add(global_id_set_->id(EntityRep<0>(index, true)),
                           ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
    }
    cell_indexset_.endResize();
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

#if HAVE_MPI

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
        if(*i<std::numeric_limits<int>::max())
            ids[*i]=idSet.id(EntityRep<codim>(i-ibegin,true));
    }
    return count;
}

template<class T>
struct GetRowType
{};

template<class T>
struct GetRowType<Opm::SparseTable<T> >
{
    typedef typename Opm::SparseTable<T>::row_type type;
};
template<class E, class A>
struct GetRowType<std::vector<E,A> >
{
    typedef typename std::vector<E,A>::value_type type;
};

PartitionType getPartitionType(const PartitionTypeIndicator& p, const EntityRep<1>& f,
                               const CpGridData&)
{
    return p.getPartitionType(f);
}

PartitionType getPartitionType(const PartitionTypeIndicator& p, int i,
                               const CpGridData& grid)
{
    return p.getPartitionType(Entity<3>(grid, i, true));
}

int getIndex(const int* i)
{
    return *i;
}

template<class T>
int getIndex(T i)
{
    return i->index();
}

struct Cell2PointsDataHandle
{
    using DataType = int;
    using Vector = std::vector<std::array<int,8> >;
    Cell2PointsDataHandle(const Vector& globalCell2Points,
                          Vector& localCell2Points,
                          std::vector<int>& flatGlobalPoints)
        : globalCell2Points_(globalCell2Points), localCell2Points_(localCell2Points),
          flatGlobalPoints_(flatGlobalPoints)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 0;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 8;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        auto i = t.index();
        const auto& points = globalCell2Points_[i];
        std::for_each(points.begin(), points.end(), [&buffer](const int& i){buffer.write(i);});
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t )
    {
        auto i = t.index();
        auto& points = localCell2Points_[i];
        std::for_each(points.begin(), points.end(),
                      [&buffer, this](int& i){
                          buffer.read(i);
                          this->flatGlobalPoints_.push_back(i);
                      });
    }
private:
    const Vector& globalCell2Points_;
    Vector& localCell2Points_;
    std::vector<int>& flatGlobalPoints_;
};

template<int from, int to>
struct RowSizeDataHandle
{
    using DataType = int;
    using Table = OrientedEntityTable<from, to>;
    RowSizeDataHandle(const Table& global,
                      std::vector<int>& noEntries)
        : global_(global), noEntries_(noEntries)
    {}
    bool fixedsize()
    {
        return true;
    }
    std::size_t size()
    {
        return 1;
    }
    template<class B>
    void gather(B& buffer, std::size_t i)
    {
        buffer.write(global_.rowSize(EntityRep<from>(i, true)));
    }
    template<class B>
    void scatter(B& buffer, std::size_t i)
    {
        buffer.read(noEntries_[i]);
    }
private:
    const Table& global_;
    std::vector<int>& noEntries_;
};

template<int from, int to>
struct OrientedEntityTableDataHandle
{
    using DataType = int;
    using Table = OrientedEntityTable<from, to>;
    OrientedEntityTableDataHandle(const Table& global,
                                  Table& local)
        : global_(global), local_(local)
    {}
    bool fixedsize()
    {
        return false;
    }
    std::size_t size(int i)
    {
        return global_.rowSize(EntityRep<from>(i, true));
    }
    template<class B>
    void gather(B& buffer, std::size_t i)
    {
        const auto& entries = global_[EntityRep<0>(i, true)];
        std::for_each(entries.begin(), entries.end(), [&buffer](const int& i){buffer.write(i);});
    }
    template<class B>
    void scatter(B& buffer, std::size_t i)
    {
        auto& entries = local_[i];
        std::for_each(entries.begin(), entries.end(), [&buffer](int& i){buffer.read(i);});
    }
private:
    const Table& global_;
    Table local_;
};

template<class T>
struct AttributeDataHandle
{
    typedef std::pair<int,char> DataType;

    AttributeDataHandle(int rank, const PartitionTypeIndicator& indicator,
                        std::vector<std::map<int, char> >& vals,
                        const T& cell_to_entity,
                        const CpGridData& grid)
        : rank_(rank), indicator_(indicator), vals_(vals),
        c2e_(cell_to_entity), grid_(grid)
    {}
    bool fixedsize()
    {
        return true;
    }
    std::size_t size(std::size_t i)
    {
        return c2e_[i].size();
    }
    template<class B>
    void gather(B& buffer, std::size_t i)
    {
        typedef typename GetRowType<T>::type::const_iterator RowIter;
        for(RowIter f=c2e_[i].begin(), fend=c2e_[i].end();
            f!=fend; ++f)
        {
            char t=getPartitionType(indicator_, *f, grid_);
            buffer.write(std::make_pair(rank_, t));
        }
    }

    template<class B>
    void scatter(B& buffer, std::size_t i, std::size_t s)
    {
        typedef typename GetRowType<T>::type::const_iterator RowIter;
        for(RowIter f=c2e_[i].begin(), fend=c2e_[i].end();
            f!=fend; ++f, --s)
        {
            std::pair<int,char> rank_attr;
            buffer.read(rank_attr);
            vals_[getIndex(f)].insert(rank_attr);
        }
    }
    int rank_;
    const PartitionTypeIndicator& indicator_;
    std::vector<std::map<int, char> >& vals_;
    const T& c2e_;
    const CpGridData& grid_;
};


template<class T, class Functor, class FromSet, class ToSet>
struct InterfaceFunctor
{
    InterfaceFunctor(std::map<int,std::pair<T,T> >& m)
        : map_(m)
    {}
    void operator()(int rank, std::size_t index, PartitionType mine, PartitionType other)
    {
        if(from.contains(mine) && to.contains(other))
            func(map_[rank].first, index);
        if(from.contains(other) && to.contains(mine))
            func(map_[rank].second, index);
    }
    std::map<int,std::pair<T,T> >& map_;
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
    void operator()(InterfaceInformation& info, std::size_t index)
    {
        info.add(index);
    }
};

template<class Tuple>
struct InterfaceTupleFunctor
{
    InterfaceTupleFunctor(Tuple& t)
    : t_(t)
    {}

    void operator()(int rank, std::size_t index, PartitionType mine, PartitionType other)
    {
        std::get<0>(t_)(rank, index, mine, other);
        std::get<1>(t_)(rank, index, mine, other);
        std::get<2>(t_)(rank, index, mine, other);
        std::get<3>(t_)(rank, index, mine, other);
        std::get<4>(t_)(rank, index, mine, other);
    }
    Tuple& t_;
};

struct Converter
{
    typedef EnumItem<PartitionType, InteriorEntity> Interior;
    typedef EnumItem<PartitionType, BorderEntity> Border;
    typedef EnumItem<PartitionType, OverlapEntity> Overlap;
    typedef EnumItem<PartitionType, FrontEntity> Front;
    typedef Combine<Interior,Border> InteriorBorder;
    typedef Combine<Overlap,Front> OverlapFront;

    typedef std::tuple<InteriorBorder, InteriorBorder,        Overlap     , Overlap,
                       AllSet<PartitionType> > SourceTuple;
    typedef std::tuple<InteriorBorder, AllSet<PartitionType>, OverlapFront, AllSet<PartitionType>,
                       AllSet<PartitionType> > DestinationTuple;
};

/**
 * \brief Reserves space for an interface.
 * \tparam i The index of the interface.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
template<std::size_t i, class InterfaceMap>
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
             interfaces)
{
    typedef typename std::map<int,std::pair<std::size_t,std::size_t> >::const_iterator Iter;
    const std::map<int,std::pair<std::size_t,std::size_t> >& sizeMap=sizes[i];
    InterfaceMap& interfaceMap=std::get<i>(interfaces);

    for(Iter iter=sizeMap.begin(), end =sizeMap.end(); iter!=end; ++iter)
    {
        std::pair<InterfaceInformation,InterfaceInformation>& interface=interfaceMap[iter->first];
        interface.first.reserve(iter->second.first);
        interface.second.reserve(iter->second.second);
    }
}

/**
 * \brief Reserves space for the interfaces.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
template<class InterfaceMap>
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
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
                            typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                            typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<std::size_t, InterfaceIncrementor,
                             typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                             typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
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
                             typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                             typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<InterfaceInformation, InterfaceAdder,
                             typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                             typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
    Base;
    AddFunctor(std::map<int,std::pair<InterfaceInformation,InterfaceInformation> >& m)
        : Base(m)
    {}
};

class FacePartitionTypeIterator
{
public:
    FacePartitionTypeIterator(const PartitionTypeIndicator* part)
    : indicator_(part), index_()
    {}
    void operator++()
    {
        ++index_;
    }
    PartitionType operator*()
    {
        return indicator_->getFacePartitionType(index_);
    }
private:
    const PartitionTypeIndicator* indicator_;
    int index_;
};


/**
 * \brief Applies a functor the each pair of the interface.
 * \tparam Functor The type of the functor to apply.
 * \param attributes[in] A vector that contains for each index a map from other
 * process ranks to the attribute there.
 * \param my_attributes[in] A vector with the attributes of each index on this process.
 * \param func The functor.
 */
template<class Functor, class T>
void iterate_over_attributes(std::vector<std::map<int,char> >& attributes,
                             T my_attribute_iter, Functor& func)
{
    typedef typename std::vector<std::map<int,char> >::const_iterator Iter;
    for(Iter begin=attributes.begin(), i=begin, end=attributes.end(); i!=end;
        ++i, ++my_attribute_iter)
    {
        typedef typename std::map<int,char>::const_iterator MIter;
        for(MIter m=i->begin(), mend=i->end(); m!=mend; ++m)
        {
            func(m->first,i-begin,PartitionType(*my_attribute_iter), PartitionType(m->second));
        }
    }
}


/**
 * \brief Creates the communication interface for either faces or points.
 * \param attributes[in] A vector that contains for each index a map from other
 * process ranks to the attribute there.
 * \param my_attributes[in] A vector with the attributes of each index on this process.
 * \param[out] interfaces The tuple with the interface maps for communication.
 */
template<class InterfaceMap,class T>
void createInterfaces(std::vector<std::map<int,char> >& attributes,
                      T partition_type_iterator,
                      std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                      interfaces)
{
    // calculate sizes
    std::vector<std::map<int,std::pair<std::size_t,std::size_t> > > sizes(5);
    typedef std::tuple<SizeFunctor<0>,SizeFunctor<1>,SizeFunctor<2>,SizeFunctor<3>,
                       SizeFunctor<4> > SizeTuple;

    SizeTuple
    size_functor_tuple = std::make_tuple(SizeFunctor<0>(sizes[0]),
                                         SizeFunctor<1>(sizes[1]),
                                         SizeFunctor<2>(sizes[2]),
                                         SizeFunctor<3>(sizes[3]),
                                         SizeFunctor<4>(sizes[4]));
    InterfaceTupleFunctor<SizeTuple> size_functor(size_functor_tuple);
    iterate_over_attributes(attributes, partition_type_iterator, size_functor);
    // reserve space
    reserve(sizes, interfaces);
    // add indices to the interface
    typedef std::tuple<AddFunctor<0>,AddFunctor<1>,AddFunctor<2>,AddFunctor<3>,
                       AddFunctor<4> > AddTuple;
    AddTuple
        add_functor_tuple(AddFunctor<0>(std::get<0>(interfaces)),
                          AddFunctor<1>(std::get<1>(interfaces)),
                          AddFunctor<2>(std::get<2>(interfaces)),
                          AddFunctor<3>(std::get<3>(interfaces)),
                          AddFunctor<4>(std::get<4>(interfaces)));
    InterfaceTupleFunctor<AddTuple> add_functor(add_functor_tuple);
    iterate_over_attributes(attributes, partition_type_iterator, add_functor);

}

#endif // #if HAVE_MPI

std::map<int,int> computeCell2Point(CpGrid& grid,
                                    const std::vector<std::array<int,8> >& globalCell2Points,
                                    std::vector<std::array<int,8> >& cell2Points,
                                    std::vector<int>& map2Global,
                                    std::size_t noCells)
{
    cell2Points.resize(noCells);
    map2Global.reserve(noCells*8);
    Cell2PointsDataHandle handle(globalCell2Points, cell2Points,
                                 map2Global);
    grid.scatterData(handle);
    // make map2Global a map from local to global
    std::sort(map2Global.begin(),map2Global.end());
    auto newEnd = std::unique(map2Global.begin(),map2Global.end());
    map2Global.resize(newEnd - map2Global.begin());
    // Convert point ids to local ones
    std::map<int, int> map2Local;
    auto current_global = map2Global.begin();
    int localId = 0;
    // \todo improve since we are inserting values by sorted keys.
    std::generate_n(std::inserter(map2Local, map2Local.begin()),
                    newEnd - current_global,
                    [&localId, &current_global](){ return std::make_pair(*(current_global++), localId++); });
    for (auto&& points : cell2Points)
    {
        for (auto&& point : points)
        {
            point = map2Local[point];
        }
    }
    return map2Local;
}

void CpGridData::distributeGlobalGrid(CpGrid& grid,
                                      const CpGridData& view_data,
                                      const std::vector<int>& cell_part)
{
#if HAVE_MPI
    // setup the remote indices.
    cell_remote_indices_.setIndexSets(cell_indexset_, cell_indexset_, ccobj_);
    cell_remote_indices_.template rebuild<false>(); // We could probably also compute this on our own, like before?

    // We can identify existing cells with the help of the index set.
    // Now we need to compute the existing faces and points. Either exist
    // if they are reachable from an existing cell.
    // We use std::numeric_limits<int>::max() to indicate non-existent entities.
    std::vector<int> map2GlobalFaceId;
    std::vector<int> map2GlobalPointId;
    std::map<int,int> point_indicator =
        computeCell2Point(grid, view_data.cell_to_point_, cell_to_point_,
                          map2GlobalPointId, cell_indexset_.size());

    map2GlobalFaceId.reserve((6*cell_indexset_.size())*1.1);

    for(ParallelIndexSet::iterator i=cell_indexset_.begin(), end=cell_indexset_.end();
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
            map2GlobalFaceId.push_back(findex);
        }
    }

    auto noExistingPoints = map2GlobalPointId.size();
    std::sort(map2GlobalFaceId.begin(),map2GlobalFaceId.end());
    auto newFaceEnd = std::unique(map2GlobalFaceId.begin(),map2GlobalFaceId.end());
    map2GlobalFaceId.resize(newFaceEnd - map2GlobalFaceId.begin());
    auto noExistingFaces = map2GlobalFaceId.size();

    // renumber  face and point indicators
    std::map<int,int> face_indicator;
    auto current_global_face = map2GlobalFaceId.begin();
    int localId = 0;
    std::generate_n(std::inserter(face_indicator, face_indicator.begin()),
                    newFaceEnd - current_global_face,
                    [&localId, &current_global_face](){ return std::make_pair(*(current_global_face++), localId++); });

    std::vector<int> map2GlobalCellId(cell_indexset_.size());
    for(ParallelIndexSet::const_iterator i=cell_indexset_.begin(), end=cell_indexset_.end();
        i!=end; ++i)
    {
        map2GlobalCellId[i->local()]=view_data.local_id_set_->id(EntityRep<0>(i->global(), true));
    }

    // Turn local face/point ids into global ones
    auto func = [&view_data](int i){
                    return view_data.local_id_set_->id(EntityRep<3>(i, true));
                };
    std::transform(map2GlobalPointId.begin(), map2GlobalPointId.end(), map2GlobalPointId.begin(),
                   func);
    auto ffunc = [&view_data](int i){
                    return view_data.local_id_set_->id(EntityRep<1>(i, true));
                };
    std::transform(map2GlobalFaceId.begin(), map2GlobalFaceId.end(), map2GlobalFaceId.begin(),
                   ffunc);
    global_id_set_->swap(map2GlobalCellId, map2GlobalFaceId, map2GlobalPointId);

    // Create the topology information. This is stored in sparse matrix like data structures.
    // First conunt the size of the nonzeros of the cell_to_face data.
    int data_size=0;
    for(auto i=cell_indexset_.begin(), end=cell_indexset_.end();
        i!=end; ++i)
    {
        const Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        data_size+=c2f.rowSize(i->global());
    }

    //- cell_to_face_ : extract owner/overlap rows from cell_to_face_
    // Construct the sparse matrix like data structure.
    //OrientedEntityTable<0, 1> cell_to_face;
    cell_to_face_.reserve(cell_indexset_.size(), data_size);

    for(auto i=cell_indexset_.begin(), end=cell_indexset_.end();
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
    for(const auto& f: face_indicator)
        data_size += f2c.rowSize(f.first);

    face_to_cell_.reserve(noExistingFaces, data_size);

    //- face_to cell_ : extract rows that connect to an existent cell
    std::vector<int> cell_indicator(view_data.cell_to_face_.size(),
                                    std::numeric_limits<int>::max());
    for(auto i=cell_indexset_.begin(), end=cell_indexset_.end(); i!=end; ++i)
        cell_indicator[i->global()]=i->local();

    for(const auto& f: face_indicator)
    {
        // face does exist
        std::vector<EntityRep<0> > new_row;
        auto old_row = f2c[f.first];
        new_row.reserve(old_row.size());
        // push back connected existent cells.
        // for those cells we use the new cell_indicator and copy the
        // orientation of the old cell.
        for(auto cell = old_row.begin(), cend=old_row.end(); cell != cend; ++cell)
        {
            // The statement below results in all faces having two neighbours
            // except for those at the domain boundary.
            // Note that along the front partition there are invalid neighbours
            // marked with index std::numeric_limits<int>::max()
            // Still they inherit the orientation to make CpGrid::faceCell happy
            new_row.push_back(EntityRep<0>(cell_indicator[cell->index()],
                                           cell->orientation()));
        }
        face_to_cell_.appendRow(new_row.begin(), new_row.end());
    }

    face_to_point_.reserve(noExistingFaces, data_size);

    //- face_to_point__ : extract row associated with existing faces_
    for(const auto& f: face_indicator)
    {
        // face does exist
        std::vector<int> new_row;
        auto old_row = view_data.face_to_point_[f.first];
        new_row.reserve(old_row.size());
        for(auto point = old_row.begin(), pend=old_row.end(); point!=pend; ++point)
        {
            new_row.push_back(point_indicator[*point]);
        }
        face_to_point_.appendRow(new_row.begin(), new_row.end());
    }

    logical_cartesian_size_=view_data.logical_cartesian_size_;

    // Set up the new topology arrays
    // Count the existing points and allocate space
    EntityVariable<cpgrid::Geometry<0, 3>, 3>& point_geom = geometry_.geomVector(std::integral_constant<int,3>());
    const std::vector<cpgrid::Geometry<0, 3> >& global_point_geom=view_data.geomVector<3>();
    point_geom.reserve(noExistingPoints);

    // Now copy the point geometries that do exist.
    for(auto pi = point_indicator.begin(), end = point_indicator.end(); pi != end;
        ++pi)
    {
        assert(pi->second == (int) point_geom.size());
        point_geom.emplace_back(global_point_geom[pi->first]);
    }


    EntityVariable<cpgrid::Geometry<3, 3>, 0>& cell_geom = geometry_.geomVector(std::integral_constant<int,0>());
    std::vector<cpgrid::Geometry<3, 3> > tmp_cell_geom(cell_indexset_.size());
    auto global_cell_geom=view_data.geomVector<0>();
    global_cell_.resize(cell_indexset_.size());
    cell_geom.resize(cell_indexset_.size());
    // Copy the existing cells.
    for (auto i = cell_indexset_.begin(), end = cell_indexset_.end(); i != end; ++i)
    {
        const auto& geom = global_cell_geom.get(i->global());
        cell_geom.get(i->local()) = Geometry<3,3>(geom.center(), geom.volume(),
                                                  point_geom,
                                                  cell_to_point_[i->local()].data());
        global_cell_[i->local()]=view_data.global_cell_[i->global()];
    }

    // count the existing faces, renumber, and allocate space.
    EntityVariable<cpgrid::Geometry<2, 3>, 1>&  face_geom = geometry_.geomVector(std::integral_constant<int,1>());
    face_geom.reserve(noExistingFaces);
    std::vector<enum face_tag> tmp_face_tag(noExistingFaces);
    std::vector<PointType> tmp_face_normals(noExistingFaces);
    const std::vector<cpgrid::Geometry<2, 3> >& global_face_geom=view_data.geomVector<1>();
    const std::vector<enum face_tag>& global_face_tag=view_data.face_tag_;
    const std::vector<PointType>& global_face_normals=view_data.face_normals_;

    // Now copy the face geometries that do exist.
    auto ft = tmp_face_tag.begin();
    auto fn = tmp_face_normals.begin();
    for (auto begin = face_indicator.begin(), fi = begin,  end = face_indicator.end(); fi != end;
        ++fi)
    {
        assert(fi->second == (int) face_geom.size());
        face_geom.push_back(global_face_geom[fi->first]);
        *ft=global_face_tag[fi->first];
        *fn=global_face_normals[fi->first];
        ++ft; ++fn;
    }
    static_cast<std::vector<PointType>&>(face_normals_).swap(tmp_face_normals);
    static_cast<std::vector<enum face_tag>&>(face_tag_).swap(tmp_face_tag);

    // - unique_boundary_ids_ : extract the ones that correspond existent faces
    if(view_data.unique_boundary_ids_.size())
    {
        // Unique boundary ids are inherited from the global grid.
        unique_boundary_ids_.reserve(view_data.face_to_cell_.size());
        for(const auto& f: face_indicator)
        {
            unique_boundary_ids_.push_back(view_data.unique_boundary_ids_[EntityRep<1>(f.first, true)]);
        }
    }
    // Compute the partition type for cell
    partition_type_indicator_->cell_indicator_.resize(cell_indexset_.size());
    for(ParallelIndexSet::const_iterator i=cell_indexset_.begin(), end=cell_indexset_.end();
            i!=end; ++i)
    {
        partition_type_indicator_->cell_indicator_[i->local()]=
            i->local().attribute()==AttributeSet::owner?
            InteriorEntity:OverlapEntity;
    }

    // Compute partition type for points
    // We initialize all points with interior. Then we loop over the faces. If a face is of
    // type border, then the type of the point is overwritten with border. In the other cases
    // we set the type of the point to the one of the face as long as the type of the point is
    // not border.
    partition_type_indicator_->point_indicator_.resize(geometry_.geomVector<3>().size(),
                                                       OverlapEntity);
    for(int i=0; i<face_to_point_.size(); ++i)
    {
        for(auto p=face_to_point_[i].begin(),
                pend=face_to_point_[i].end(); p!=pend; ++p)
        {
            PartitionType new_type=partition_type_indicator_->getFacePartitionType(i);
            PartitionType old_type=PartitionType(partition_type_indicator_->point_indicator_[*p]);
            if(old_type==InteriorEntity)
            {
                if(new_type!=OverlapEntity)
                    partition_type_indicator_->point_indicator_[*p]=new_type;
            }
            if(old_type==OverlapEntity)
                partition_type_indicator_->point_indicator_[*p]=new_type;
            if(old_type==FrontEntity && new_type==BorderEntity)
                partition_type_indicator_->point_indicator_[*p]=new_type;
        }
    }

    // Compute the interface information for cells
    std::get<InteriorBorder_All_Interface>(cell_interfaces_)
        .build(cell_remote_indices_, EnumItem<AttributeSet, AttributeSet::owner>(),
               AllSet<AttributeSet>());
    std::get<Overlap_OverlapFront_Interface>(cell_interfaces_)
        .build(cell_remote_indices_, EnumItem<AttributeSet, AttributeSet::copy>(),
               EnumItem<AttributeSet, AttributeSet::copy>());
    std::get<Overlap_All_Interface>(cell_interfaces_)
        .build(cell_remote_indices_, EnumItem<AttributeSet, AttributeSet::copy>(),
                                 AllSet<AttributeSet>());
    std::get<All_All_Interface>(cell_interfaces_)
        .build(cell_remote_indices_, AllSet<AttributeSet>(), AllSet<AttributeSet>());

    // Now we use the all_all communication of the cells to compute which faces and points
    // are also present on other processes and with what attribute.
    const auto& all_all_cell_interface = std::get<All_All_Interface>(cell_interfaces_);

    // Work around a bug/deadlock in DUNE <=2.5.1 which happens if the
    // buffer cannot hold all data that needs to be send.
    // https://gitlab.dune-project.org/core/dune-common/merge_requests/416
    // For this we calculate an upper barrier of the number of
    // data items to be send manually and use it to construct a
    // VariableSizeCommunicator with sufficient buffer.
    std::size_t max_entries = 0;
    for (const auto& pair: all_all_cell_interface.interfaces() )
    {
        using std::max;
        max_entries = max(max_entries, pair.second.first.size());
        max_entries = max(max_entries, pair.second.second.size());
    }
    Dune::VariableSizeCommunicator<> comm(all_all_cell_interface.communicator(),
                                          all_all_cell_interface.interfaces(),
                                          max_entries*8*sizeof(int));
    /*
      // code deactivated, because users cannot access face indices and therefore
      // communication on faces makes no sense!
    std::vector<std::map<int,char> > face_attributes(noExistingFaces);
    AttributeDataHandle<Opm::SparseTable<EntityRep<1> > >
        face_handle(ccobj_.rank(), *partition_type_indicator_,
                    face_attributes, static_cast<Opm::SparseTable<EntityRep<1> >&>(cell_to_face_),
                    *this);
    if( std::get<All_All_Interface>(cell_interfaces_).interfaces().size() )
    {
        comm.forward(face_handle);
    }
    createInterfaces(face_attributes, FacePartitionTypeIterator(partition_type_indicator_),
                     face_interfaces_);
    std::vector<std::map<int,char> >().swap(face_attributes);
    */
    std::vector<std::map<int,char> > point_attributes(noExistingPoints);
    AttributeDataHandle<std::vector<std::array<int,8> > >
        point_handle(ccobj_.rank(), *partition_type_indicator_,
                     point_attributes, cell_to_point_, *this);
    if( static_cast<const Dune::Interface&>(std::get<All_All_Interface>(cell_interfaces_))
        .interfaces().size() )
    {
        comm.forward(point_handle);
    }
    createInterfaces(point_attributes, partition_type_indicator_->point_indicator_.begin(),
                     point_interfaces_);
#else // #if HAVE_MPI
    static_cast<void>(grid);
    static_cast<void>(view_data);
    static_cast<void>(cell_part);
    static_cast<void>(overlap_layers);
#endif
}

} // end namespace cpgrid
} // end namespace Dune
