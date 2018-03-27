#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H

#include <cstdint>
#include <vector>


namespace PHARE
{
//! base class for NdArrayVector 1D, 2D and 3D.
/**
 * This base class gathers all code that is common to 1D, 2D and 3D implementations.
 */
template<typename DataType = double>
class NdArrayVectorBase
{
protected:
    NdArrayVectorBase() = delete;
    explicit NdArrayVectorBase(std::size_t size)
        : data_(size)
    {
    }

    NdArrayVectorBase(NdArrayVectorBase const& source) = default;
    NdArrayVectorBase(NdArrayVectorBase&& source)      = default;
    NdArrayVectorBase& operator=(NdArrayVectorBase const& source) = default;
    NdArrayVectorBase& operator=(NdArrayVectorBase&& source) = default;

    std::vector<DataType> data_;

public:
    //! user can check data_type to know of which type the elements are
    using data_type = DataType;

    //! return the total number of elements in the container
    std::size_t size()
    {
        auto s = data_.size();
        return static_cast<std::size_t>(s);
    }
};



//! NdArrayVector1D is a 1D data container implementation
/** NdArrayVector1D uses a contiguous std::vector as its internal
 *  data representation. It can store any kind of elements although 'double'
 *  is the default type.
 */
template<typename DataType = double>
class NdArrayVector1D : public NdArrayVectorBase<DataType>
{
public:
    //! builds an NdArrayVector1D by specifying its number of elements
    explicit NdArrayVector1D(int nx)
        : NdArrayVectorBase<DataType>(nx)
        , nx_{nx}
    {
    }

    NdArrayVector1D()                              = delete;
    NdArrayVector1D(NdArrayVector1D const& source) = default;
    NdArrayVector1D(NdArrayVector1D&& source)      = default;
    NdArrayVector1D& operator=(NdArrayVector1D const& source) = default;
    NdArrayVector1D& operator=(NdArrayVector1D&& source) = default;

    //! read/write access operator
    DataType& operator()(int i) { return this->data_[i]; }
    DataType const& operator()(int i) const { return this->data_[i]; }

    static const int dimension = 1;

private:
    int nx_ = 0;
};


//! NdArrayVector2D is an implementation for a 2-dimensional container
/** This representation offers a contiguous representation which can be
 *  efficient for repetitive random access to the memory.
 *  Elements are stored following the C order, i.e. in array(i,j), 'i' is the
 *  slower varying index.
 */
template<typename DataType = double>
class NdArrayVector2D : public NdArrayVectorBase<DataType>
{
public:
    //! build a NdArrayVector2D specifying its number of elements in the 1st and 2nd dims.
    NdArrayVector2D(int nx, int ny)
        : NdArrayVectorBase<DataType>(nx * ny)
        , nx_{nx}
        , ny_{ny}
    {
    }

    NdArrayVector2D()                              = delete;
    NdArrayVector2D(NdArrayVector2D const& source) = default;
    NdArrayVector2D(NdArrayVector2D&& source)      = default;
    NdArrayVector2D& operator=(NdArrayVector2D const& source) = default;
    NdArrayVector2D& operator=(NdArrayVector2D&& source) = default;

    //! read/write data access operator returns C-ordered data.
    DataType& operator()(int i, int j) { return this->data_[j + ny_ * i]; }

    //! read only access. See read/write.
    DataType const& operator()(int i, int j) const { return this->data_[linearIt(i, j)]; }

    static const int dimension = 2;

private:
    int constexpr linearIt(int i, int j) const { return j + ny_ * i; }


    int nx_ = 0;
    int ny_ = 0;
};



//! NdArrayVector3D is an implementation for a 3-dimensional container
/** behaves as the 2D version.
 */
template<typename DataType = double>
class NdArrayVector3D : public NdArrayVectorBase<DataType>
{
public:
    NdArrayVector3D(int nx, int ny, int nz)
        : NdArrayVectorBase<DataType>(nx * ny * nz)
        , nx_{nx}
        , ny_{ny}
        , nz_{nz}
    {
    }


    NdArrayVector3D()                              = delete;
    NdArrayVector3D(NdArrayVector3D const& source) = default;
    NdArrayVector3D(NdArrayVector3D&& source)      = default;
    NdArrayVector3D& operator=(NdArrayVector3D const& source) = default;
    NdArrayVector3D& operator=(NdArrayVector3D&& source) = default;


    DataType& operator()(int i, int j, int k) { return this->data_[k + nz_ * (j + ny_ * i)]; }
    DataType const& operator()(int i, int j, int k) const { return this->data_[linearIt(i, j, k)]; }

    static const int dimension = 3;

private:
    int constexpr linearIt(int i, int j, int k) const { return k + nz_ * (j + ny_ * i); }
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
};


} // namespace PHARE
#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
