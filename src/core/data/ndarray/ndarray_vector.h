#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H

#include <array>
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



    auto begin() const { return std::begin(data_); }
    auto begin() { return std::begin(data_); }

    auto end() { return std::end(data_); }
    auto end() const { return std::end(data_); }


    void zero()
    {
        for (auto& v : data_)
        {
            v = DataType{0};
        }
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
    explicit NdArrayVector1D(uint32_t nx)
        : NdArrayVectorBase<DataType>(nx)
        , nx_{nx}
    {
    }

    explicit NdArrayVector1D(std::array<uint32_t, 1> const& nCell)
        : NdArrayVectorBase<DataType>(nCell[0])
        , nx_{nCell[0]}
    {
    }

    NdArrayVector1D()                              = delete;
    NdArrayVector1D(NdArrayVector1D const& source) = default;
    NdArrayVector1D(NdArrayVector1D&& source)      = default;
    NdArrayVector1D& operator                      =(NdArrayVector1D const& source)
    {
        if (nx_ != source.nx_)
        {
            throw std::runtime_error(
                "Error NdArrayVector1D cannot be assigned, incompatible sizes");
        }
        else
        {
            this->data_ = source.data_;
        }
        return *this;
    }

    NdArrayVector1D& operator=(NdArrayVector1D&& source)
    {
        if (nx_ != source.nx_)
        {
            throw std::runtime_error(
                "Error NdArrayVector1D cannot be assigned, incompatible sizes");
        }
        else
        {
            this->data_ = std::move(source.data_);
        }
        return *this;
    }

    //! read/write access operator
    DataType& operator()(uint32_t i) { return this->data_[i]; }
    DataType const& operator()(uint32_t i) const { return this->data_[i]; }

    static const int dimension = 1;
    using type                 = DataType;

private:
    uint32_t nx_ = 0;
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
    NdArrayVector2D(uint32_t nx, uint32_t ny)
        : NdArrayVectorBase<DataType>(nx * ny)
        , nx_{nx}
        , ny_{ny}
    {
    }

    explicit NdArrayVector2D(std::array<uint32_t, 2> const& nbCell)
        : NdArrayVectorBase<DataType>(nbCell[0] * nbCell[1])
        , nx_{nbCell[0]}
        , ny_{nbCell[1]}
    {
    }

    NdArrayVector2D()                              = delete;
    NdArrayVector2D(NdArrayVector2D const& source) = default;
    NdArrayVector2D(NdArrayVector2D&& source)      = default;
    NdArrayVector2D& operator                      =(NdArrayVector2D const& source)
    {
        if (nx_ != source.nx_ || ny_ != source.ny_)
        {
            throw std::runtime_error(
                "Error NdArrayVector1D cannot be assigned, incompatible sizes");
        }
        else
        {
            this->data_ = source.data_;
        }
        return *this;
    }

    NdArrayVector2D& operator=(NdArrayVector2D&& source)
    {
        if (nx_ != source.nx_ || ny_ != source.ny_)
        {
            throw std::runtime_error(
                "Error NdArrayVector1D cannot be assigned, incompatible sizes");
        }
        else
        {
            this->data_ = std::move(source.data_);
        }
        return *this;
    }

    //! read/write data access operator returns C-ordered data.
    DataType& operator()(uint32_t i, uint32_t j) { return this->data_[j + ny_ * i]; }

    //! read only access. See read/write.
    DataType const& operator()(uint32_t i, uint32_t j) const { return this->data_[linearIt(i, j)]; }

    static const int dimension = 2;
    using type                 = DataType;

private:
    int constexpr linearIt(int i, int j) const { return j + ny_ * i; }


    uint32_t nx_ = 0;
    uint32_t ny_ = 0;
};



//! NdArrayVector3D is an implementation for a 3-dimensional container
/** behaves as the 2D version.
 */
template<typename DataType = double>
class NdArrayVector3D : public NdArrayVectorBase<DataType>
{
public:
    NdArrayVector3D(uint32_t nx, uint32_t ny, uint32_t nz)
        : NdArrayVectorBase<DataType>(nx * ny * nz)
        , nx_{nx}
        , ny_{ny}
        , nz_{nz}
    {
    }

    explicit NdArrayVector3D(std::array<uint32_t, 3> const& nbCell)
        : NdArrayVectorBase<DataType>(nbCell[0] * nbCell[1] * nbCell[2])
        , nx_{nbCell[0]}
        , ny_{nbCell[1]}
        , nz_{nbCell[2]}
    {
    }


    NdArrayVector3D()                              = delete;
    NdArrayVector3D(NdArrayVector3D const& source) = default;
    NdArrayVector3D(NdArrayVector3D&& source)      = default;
    NdArrayVector3D& operator                      =(NdArrayVector3D const& source)
    {
        if (nx_ != source.nx_ || ny_ != source.ny_ || nz_ != source.nz_)
        {
            throw std::runtime_error(
                "Error NdArrayVector1D cannot be assigned, incompatible sizes");
        }
        else
        {
            this->data_ = source.data_;
        }
        return *this;
    }

    NdArrayVector3D& operator=(NdArrayVector3D&& source)
    {
        if (nx_ != source.nx_ || ny_ != source.ny_ || nz_ != source.nz_)
        {
            throw std::runtime_error(
                "Error NdArrayVector1D cannot be assigned, incompatible sizes");
        }
        else
        {
            this->data_ = std::move(source.data_);
        }
        return *this;
    }


    DataType& operator()(uint32_t i, uint32_t j, uint32_t k)
    {
        return this->data_[k + nz_ * (j + ny_ * i)];
    }
    DataType const& operator()(uint32_t i, uint32_t j, uint32_t k) const
    {
        return this->data_[linearIt(i, j, k)];
    }

    static const int dimension = 3;
    using type                 = DataType;

private:
    int constexpr linearIt(uint32_t i, uint32_t j, uint32_t k) const
    {
        return k + nz_ * (j + ny_ * i);
    }
    uint32_t nx_ = 0;
    uint32_t ny_ = 0;
    uint32_t nz_ = 0;
};


} // namespace PHARE
#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
