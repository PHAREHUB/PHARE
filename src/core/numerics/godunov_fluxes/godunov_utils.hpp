#ifndef CORE_NUMERICS_GODUNOV_GODUNOV_UTILS_HPP
#define CORE_NUMERICS_GODUNOV_GODUNOV_UTILS_HPP

#include "core/data/field/field.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"
#include "core/utilities/index/index.hpp"
#include <cassert>

namespace PHARE::core
{
template<typename Float = double>
struct PerIndexVector
{
    // macOS compiler does not support template deduction for aggregates
    PerIndexVector(Float x, Float y, Float z)
        : x{x}
        , y{y}
        , z{z}
    {
    }

    Float x, y, z;
};

template<typename Float>
struct PerIndex
{
    PerIndex(Float rho, PerIndexVector<Float> V, PerIndexVector<Float> B, Float P)
        : rho{rho}
        , V{V}
        , B{B}
        , P{P}
    {
    }

    auto as_tuple() const { return std::make_tuple(rho, V.x, V.y, V.z, B.x, B.y, B.z, P); }

    void to_conservative(auto const& gamma)
    {
        auto const [rhoVx, rhoVy, rhoVz] = vToRhoV(rho, V.x, V.y, V.z);
        Float Etot                       = eosPToEtot(gamma, rho, V.x, V.y, V.z, B.x, B.y, B.z, P);

        V.x = rhoVx;
        V.y = rhoVy;
        V.z = rhoVz;
        P   = Etot;
    }

    template<typename T>
    PerIndex& operator=(PerIndex<T> const& other)
    {
        rho = other.rho;
        V.x = other.V.x;
        V.y = other.V.y;
        V.z = other.V.z;
        B.x = other.B.x;
        B.y = other.B.y;
        B.z = other.B.z;
        P   = other.P;
        return *this;
    }

    auto& rhoV() { return V; }

    auto& rhoV() const { return V; }

    auto& Etot() { return P; }

    auto& Etot() const { return P; }

    Float rho;
    PerIndexVector<Float> V;
    PerIndexVector<Float> B;
    Float P;

#ifndef NDEBUG
    bool isConservative{
        true}; // does nothing, we need a better system if we want to enforce this (because the we
               // also create already conservative versions of this structure)
#endif
};

template<typename T>
struct is_field_or_tensorfield : std::false_type
{
};

template<std::size_t dim, typename PhysicalQuantity>
struct is_field_or_tensorfield<Field<dim, PhysicalQuantity>> : std::true_type
{
};

template<typename Field_t, typename PhysicalQuantity>
struct is_field_or_tensorfield<TensorField<Field_t, PhysicalQuantity>> : std::true_type
{
};

template<typename T>
concept NotFieldOrTensorField = !is_field_or_tensorfield<T>::value;

template<typename T>
concept ViewVector = requires(T t) {
    typename T::value_type;
    t[0];
} && NotFieldOrTensorField<T>;

template<typename T>
struct get_value_type
{
    using type = T;
};

template<ViewVector T>
struct get_value_type<T>
{
    using type = typename T::value_type;
};

template<typename T>
using value_type_t = typename get_value_type<T>::type;

template<typename Field, typename VecField>
struct AllFluxes
{
    using Float = typename Field::value_type;

    AllFluxes<value_type_t<Field>, value_type_t<VecField>> operator[](std::size_t i)
        requires((ViewVector<Field>) && (ViewVector<VecField>))
    {
        return AllFluxes<value_type_t<Field>, value_type_t<VecField>>{
            rho_fx[i], rhoV_fx[i], B_fx[i],   Etot_fx[i], rho_fy[i], rhoV_fy[i],
            B_fy[i],   Etot_fy[i], rho_fz[i], rhoV_fz[i], B_fz[i],   Etot_fz[i]};
    }

    template<auto direction>
    auto get_dir(MeshIndex<Field::dimension> index) const
        requires((!ViewVector<Field>) && (!ViewVector<VecField>))
    {
        using Float = typename Field::value_type;

        if constexpr (direction == Direction::X)
            return PerIndex<Float&>{
                rho_fx(index),
                {rhoV_fx(Component::X)(index), rhoV_fx(Component::Y)(index),
                 rhoV_fx(Component::Z)(index)},
                {B_fx(Component::X)(index), B_fx(Component::Y)(index), B_fx(Component::Z)(index)},
                Etot_fx(index)};
        else if constexpr (direction == Direction::Y)
            return PerIndex<Float&>{
                rho_fy(index),
                {rhoV_fy(Component::X)(index), rhoV_fy(Component::Y)(index),
                 rhoV_fy(Component::Z)(index)},
                {B_fy(Component::X)(index), B_fy(Component::Y)(index), B_fy(Component::Z)(index)},
                Etot_fy(index)};
        else if constexpr (direction == Direction::Z)
            return PerIndex<Float&>{
                rho_fz(index),
                {rhoV_fz(Component::X)(index), rhoV_fz(Component::Y)(index),
                 rhoV_fz(Component::Z)(index)},
                {B_fz(Component::X)(index), B_fz(Component::Y)(index), B_fz(Component::Z)(index)},
                Etot_fz(index)};
    }

    template<auto direction>
    auto& expose_centering() const
        requires((!ViewVector<Field>) && (!ViewVector<VecField>))
    {
        if constexpr (direction == Direction::X)
            return rho_fx;
        else if constexpr (direction == Direction::Y)
            return rho_fy;
        else if constexpr (direction == Direction::Z)
            return rho_fz;
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(rho_fx, rhoV_fx, B_fx, Etot_fx, rho_fy, rhoV_fy, B_fy, Etot_fy,
                                     rho_fz, rhoV_fz, B_fz, Etot_fz);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(rho_fx, rhoV_fx, B_fx, Etot_fx, rho_fy, rhoV_fy, B_fy, Etot_fy,
                                     rho_fz, rhoV_fz, B_fz, Etot_fz);
    }

    Field rho_fx;
    VecField rhoV_fx;
    VecField B_fx;
    Field Etot_fx;

    Field rho_fy;
    VecField rhoV_fy;
    VecField B_fy;
    Field Etot_fy;

    Field rho_fz;
    VecField rhoV_fz;
    VecField B_fz;
    Field Etot_fz;
};

} // namespace PHARE::core

#endif
