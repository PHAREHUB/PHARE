#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H

#include <cstddef>
#include <vector>

#include "particle.h"

namespace PHARE::core
{
template<std::size_t dim, bool contiguous = false>
struct ParticleArray;

template<std::size_t dim>
struct ParticleArray<dim, false>
{
    using value_type                   = Particle<dim, false>;
    using iterator                     = typename std::vector<value_type>::iterator;
    static constexpr bool is_contigous = false;

    ParticleArray() {}

    ParticleArray(size_t size)
        : particles{size}
    {
    }

    auto& operator[](size_t idx) { return particles[idx]; }
    auto& operator[](size_t idx) const { return particles[idx]; }

    template<typename T>
    auto& emplace_back()
    {
        return particles.emplace_back();
    }

    template<typename... Ts>
    auto& emplace_back(Ts&&... args)
    {
        return particles.emplace_back(std::forward<Ts>(args)...);
    }

    template<typename T>
    void push_back(T t)
    {
        particles.push_back(t);
    }

    auto begin() const { return particles.begin(); }
    auto begin() { return particles.begin(); }
    auto end() const { return particles.end(); }
    auto end() { return particles.end(); }

    void clear() { particles.clear(); }

    template<typename A, typename B>
    void erase(A a, B b)
    {
        particles.erase(a, b);
    }

    template<typename A, typename B, typename C>
    void insert(A a, B b, C c)
    {
        particles.insert(a, b, c);
    }

    auto size() const { return particles.size(); }


    std::vector<value_type> particles;

    ParticleArray(const ParticleArray&)             = delete;
    ParticleArray(const ParticleArray&&)            = delete;
    ParticleArray& operator&(const ParticleArray&)  = delete;
    ParticleArray& operator&(const ParticleArray&&) = delete;
};

template<typename type, size_t dim>
struct ParticleArrayWrapper
{
    ParticleArrayWrapper(type* wrapping)
        : pointer{wrapping}
    {
    }

    ParticleArrayWrapper& operator=(std::array<type, dim> const&& input)
    {
        return operator=(input);
    }
    ParticleArrayWrapper& operator=(std::array<type, dim> const& input)
    {
        std::copy(&input[0], &input[0] + dim, pointer);
        return *this;
    }

    type& operator[](size_t idx) const { return pointer[idx]; }
    type& operator[](size_t idx) { return pointer[idx]; }

    /*
        operator type*() { return pointer; }
        operator std::array<type, dim>() const
        {
            std::array<type, dim> array;
            std::copy(pointer, pointer + dim, array);
            return array;
        }*/

    type* pointer;
};

template<std::size_t dim>
struct ParticleArrayIndex
{
    double &weight, &charge;
    ParticleArrayWrapper<int, dim> iCell;
    ParticleArrayWrapper<float, dim> delta;
    ParticleArrayWrapper<double, 3> v;
};

template<std::size_t dim>
struct ParticleArray<dim, true> : public Particle<dim, true>
{
    using Particle<dim, true>::idx_;
    using Particle<dim, true>::size_;
    using Particle<dim, true>::weight;
    using Particle<dim, true>::charge;
    using Particle<dim, true>::iCell;
    using Particle<dim, true>::delta;
    using Particle<dim, true>::v;
    static constexpr bool is_contigous = true;

    ParticleArray() {}

    ParticleArray(size_t size)
        : Particle<dim, true>{size}
    {
    }

    auto vectorTuple() { return std::forward_as_tuple(iCell, delta, weight, charge, v); }

    auto& emplace_back()
    {
        if (idx_ == size_)
        {
            for (size_t i = 0; i < dim; i++)
            {
                iCell.emplace_back();
                delta.emplace_back();
            }
            for (size_t i = 0; i < 3; i++)
                v.emplace_back();

            charge.emplace_back();
            weight.emplace_back();
            size_++;
        }
        idx_ptr = get_array_index(idx_++);
        return *idx_ptr.get();
    }

    auto get_array_index(size_t idx)
    {
        return std::move(std::make_shared<aggregate_adapter<ParticleArrayIndex<dim>>>(
            weight[idx], charge[idx], &iCell[idx * dim], &delta[idx_ * dim], &v[idx * 3]));
    }

    class iterator
    {
    public:
        iterator(ParticleArray<dim, true>& outer, size_t _position)
            : idx_(_position)
            , outer_(outer)
        {
            set();
        }
        void set() { ptr = outer_.get_array_index(idx_); }
        auto& operator*() { return *ptr.get(); }
        auto& operator*() const { return *ptr.get(); }
        iterator& operator++()
        {
            ++idx_;
            set();
            return *this;
        }
        bool operator!=(const iterator& it) const { return idx_ != it.idx_; }

    private:
        size_t idx_;
        ParticleArray<dim, true>& outer_;
        std::shared_ptr<ParticleArrayIndex<dim>> ptr;
    };

    auto begin() const { return iterator(*this, 0); }
    auto begin() { return iterator(*this, 0); }
    auto end() const { return iterator(*this, size_); }
    auto end() { return iterator(*this, size_); }

    void clear()
    {
        core::apply(vectorTuple(), [](auto& arg) { arg.clear(); });
        size_ = idx_ = 0;
    }

    static void swap(ParticleArray& left, ParticleArray& right)
    {
        size_t idx     = 0;
        auto thisTuple = left.vectorTuple();
        core::apply(right.vectorTuple(), [&](auto& arg) {
            std::swap(std::get<idx>(thisTuple), arg);
            idx++;
        });
    }

    std::shared_ptr<ParticleArrayIndex<dim>> idx_ptr;

    ParticleArray(const ParticleArray&)             = delete;
    ParticleArray(const ParticleArray&&)            = delete;
    ParticleArray& operator&(const ParticleArray&)  = delete;
    ParticleArray& operator&(const ParticleArray&&) = delete;
};

template<std::size_t dim, bool contiguous = false>
void empty(ParticleArray<dim, contiguous>& array)
{
    array.clear();
}

template<size_t dim, bool contiguous>
void swap(ParticleArray<dim, contiguous>& array1, ParticleArray<dim, contiguous>& array2)
{
    if constexpr (contiguous)
        ParticleArray<dim, contiguous>::swap(array1, array2);
    else
        std::swap(array1.particles, array2.particles);
}

} // namespace PHARE::core

#endif
