#ifndef OPM_ENTITY2INDEXDATAHANDLE_HEADER
#define OPM_ENTITY2INDEXDATAHANDLE_HEADER
namespace Dune
{
namespace cpgrid
{
template<int codim> class Entity;

/// \brief Wrapper that turns a data handle suitable for dune-grid into one based on
/// integers instead of entities.
///
/// \tparam DataHandle The type of the data handle to wrap. Has to adhere to the interface
/// of Dune::DataHandleIf
///  \tparam codim The codimension to use when mapping indices to Entities.
template<class DataHandle, int codim>
class Entity2IndexDataHandle
{
public:
    typedef typename DataHandle::DataType DataType;
    
    Entity2IndexDataHandle(const CpGridData& grid, DataHandle& data)
        : grid_(grid), data_(data)
    {}
    bool fixedsize()
    {
        return data_.fixedsize(3, codim);
    }
    std::size_t size(std::size_t i)
    {
        return data_.size(Entity<codim>(grid_, i, true));
    }
    template<class B>
    void gather(B& buffer, std::size_t i)
    {
        data_.gather(buffer, Entity<codim>(grid_, i, true));
    }
    template<class B>
    void scatter(B& buffer, std::size_t i, std::size_t s)
    {
        data_.scatter(buffer, Entity<codim>(grid_, i, true), s);
    }
private:
    const CpGridData& grid_;
    DataHandle& data_;
    
};
} // end namespace cpgrid
} // end namespace Dune
#endif
