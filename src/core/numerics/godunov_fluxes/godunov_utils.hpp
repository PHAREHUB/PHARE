#ifndef CORE_NUMERICS_GODUNOV_GODUNOV_UTILS_HPP
#define CORE_NUMERICS_GODUNOV_GODUNOV_UTILS_HPP

#include "core/data/field/field.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/vecfield/vecfield.hpp"
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
    using Float                     = typename Field::value_type;
    static constexpr auto dimension = Field::dimension;

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
        if constexpr (dimension == 1)
            return std::forward_as_tuple(rho_fx, rhoV_fx, B_fx, Etot_fx);
        else if constexpr (dimension == 2)
            return std::forward_as_tuple(rho_fx, rhoV_fx, B_fx, Etot_fx, rho_fy, rhoV_fy, B_fy,
                                         Etot_fy);
        else if constexpr (dimension == 3)
            return std::forward_as_tuple(rho_fx, rhoV_fx, B_fx, Etot_fx, rho_fy, rhoV_fy, B_fy,
                                         Etot_fy, rho_fz, rhoV_fz, B_fz, Etot_fz);
        else
            throw std::runtime_error("Error - AllFluxes - dimension not supported");
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        if constexpr (dimension == 1)
            return std::forward_as_tuple(rho_fx, rhoV_fx, B_fx, Etot_fx);
        else if constexpr (dimension == 2)
            return std::forward_as_tuple(rho_fx, rhoV_fx, B_fx, Etot_fx, rho_fy, rhoV_fy, B_fy,
                                         Etot_fy);
        else if constexpr (dimension == 3)
            return std::forward_as_tuple(rho_fx, rhoV_fx, B_fx, Etot_fx, rho_fy, rhoV_fy, B_fy,
                                         Etot_fy, rho_fz, rhoV_fz, B_fz, Etot_fz);
        else
            throw std::runtime_error("Error - AllFluxes - dimension not supported");
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

struct AllFluxesNames
{
    std::string rho_fx;
    std::string rhoV_fx;
    std::string B_fx;
    std::string Etot_fx;

    std::string rho_fy;
    std::string rhoV_fy;
    std::string B_fy;
    std::string Etot_fy;

    std::string rho_fz;
    std::string rhoV_fz;
    std::string B_fz;
    std::string Etot_fz;

    AllFluxesNames() = default;

    template<typename AllFluxesT>
    explicit AllFluxesNames(AllFluxesT const& f)
        : rho_fx{f.rho_fx.name()}
        , rhoV_fx{f.rhoV_fx.name()}
        , B_fx{f.B_fx.name()}
        , Etot_fx{f.Etot_fx.name()}
        , rho_fy{f.rho_fy.name()}
        , rhoV_fy{f.rhoV_fy.name()}
        , B_fy{f.B_fy.name()}
        , Etot_fy{f.Etot_fy.name()}
        , rho_fz{f.rho_fz.name()}
        , rhoV_fz{f.rhoV_fz.name()}
        , B_fz{f.B_fz.name()}
        , Etot_fz{f.Etot_fz.name()}
    {
    }
};

// maybe we could want something more general than this, and use class iterators instead.
// if not, we could consider using concepts to make sure this is not used in the wrong context
template<typename Layout, typename Fn, typename First, typename... Fluxes>
void evalFluxesOnGhostBox(Layout& layout, Fn&& fn, First& first, Fluxes&... fluxes)
{
    auto static constexpr dimension = std::decay_t<decltype(layout)>::dimension;

    auto evalField = [&](auto& firstField, auto&... fluxFields) {
        layout.evalOnGhostBox(firstField, [&](auto const&... args) mutable {
            if constexpr (sizeof...(Fluxes) > 0)
            {
                fn(firstField, fluxFields..., args...);
            }
            else
            {
                fn(firstField, args...);
            }
        });
    };

    evalField(first.rho_fx, fluxes.rho_fx...);
    evalField(first.rhoV_fx(core::Component::X), fluxes.rhoV_fx(core::Component::X)...);
    evalField(first.rhoV_fx(core::Component::Y), fluxes.rhoV_fx(core::Component::Y)...);
    evalField(first.rhoV_fx(core::Component::Z), fluxes.rhoV_fx(core::Component::Z)...);

    evalField(first.B_fx(core::Component::X), fluxes.B_fx(core::Component::X)...);
    evalField(first.B_fx(core::Component::Y), fluxes.B_fx(core::Component::Y)...);
    evalField(first.B_fx(core::Component::Z), fluxes.B_fx(core::Component::Z)...);

    evalField(first.Etot_fx, fluxes.Etot_fx...);

    if constexpr (dimension >= 2)
    {
        evalField(first.rho_fy, fluxes.rho_fy...);
        evalField(first.rhoV_fy(core::Component::X), fluxes.rhoV_fy(core::Component::X)...);
        evalField(first.rhoV_fy(core::Component::Y), fluxes.rhoV_fy(core::Component::Y)...);
        evalField(first.rhoV_fy(core::Component::Z), fluxes.rhoV_fy(core::Component::Z)...);

        evalField(first.B_fy(core::Component::X), fluxes.B_fy(core::Component::X)...);
        evalField(first.B_fy(core::Component::Y), fluxes.B_fy(core::Component::Y)...);
        evalField(first.B_fy(core::Component::Z), fluxes.B_fy(core::Component::Z)...);

        evalField(first.Etot_fy, fluxes.Etot_fy...);

        if constexpr (dimension == 3)
        {
            evalField(first.rho_fz, fluxes.rho_fz...);
            evalField(first.rhoV_fz(core::Component::X), fluxes.rhoV_fz(core::Component::X)...);
            evalField(first.rhoV_fz(core::Component::Y), fluxes.rhoV_fz(core::Component::Y)...);
            evalField(first.rhoV_fz(core::Component::Z), fluxes.rhoV_fz(core::Component::Z)...);

            evalField(first.B_fz(core::Component::X), fluxes.B_fz(core::Component::X)...);
            evalField(first.B_fz(core::Component::Y), fluxes.B_fz(core::Component::Y)...);
            evalField(first.B_fz(core::Component::Z), fluxes.B_fz(core::Component::Z)...);

            evalField(first.Etot_fz, fluxes.Etot_fz...);
        }
    }
}

} // namespace PHARE::core

#endif
