#ifndef PHARE_UPWIND_CONSTRAINED_TRANSPORT_UTILS_HPP
#define PHARE_UPWIND_CONSTRAINED_TRANSPORT_UTILS_HPP

#include "core/def.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include <vector>

namespace PHARE::core
{

// Bundles heterogeneous resources (Fields and VecFields), so they can be exposed by the runtime
// resource view list of UpwindConstrainedTransportState. Hence, related memory can be allocated
// only when necessary based on a runtime criterion.
template<typename VecField>
struct UpwindConstrainedTransportTransverseState
{
    using Field                     = VecField::field_type;
    constexpr static auto dimension = VecField::dimension;

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        if constexpr (dimension == 1)
            return std::forward_as_tuple(jt_x, rhot_x);
        else if constexpr (dimension == 2)
            return std::forward_as_tuple(jt_x, rhot_x, jt_y, rhot_y);
        else
            return std::forward_as_tuple(jt_x, rhot_x, jt_y, rhot_y, jt_z, rhot_z);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        if constexpr (dimension == 1)
            return std::forward_as_tuple(jt_x, rhot_x);
        else if constexpr (dimension == 2)
            return std::forward_as_tuple(jt_x, rhot_x, jt_y, rhot_y);
        else
            return std::forward_as_tuple(jt_x, rhot_x, jt_y, rhot_y, jt_z, rhot_z);
    }

    template<typename View>
    auto as(auto&& a, auto&&... args) const
    {
        return View{{a(jt_x, args...)},   {a(jt_y, args...)},   {a(jt_z, args...)},
                    {a(rhot_x, args...)}, {a(rhot_y, args...)}, {a(rhot_z, args...)}};
    }

    /// those 'transverse' current and density fields are an average value on faces (hence _x,
    /// _y, _z) of left and right states. The averaging formula depends on the Riemann solver.
    VecField jt_x{"j_t_x", MHDQuantity::Vector::VecFlux_x};
    VecField jt_y{"j_t_y", MHDQuantity::Vector::VecFlux_y};
    VecField jt_z{"j_t_z", MHDQuantity::Vector::VecFlux_z};

    Field rhot_x{"rho_t_x", MHDQuantity::Scalar::ScalarFlux_x};
    Field rhot_y{"rho_t_y", MHDQuantity::Scalar::ScalarFlux_y};
    Field rhot_z{"rho_t_z", MHDQuantity::Scalar::ScalarFlux_z};
};


// the fields and accessors shared by UpwindConstrainedTransportState and its view
//  basic::UpwindConstrainedTransportState, Derived::transverse_state() gives the transverse state
template<typename Derived, typename VecField>
struct UpwindConstrainedTransportFields
{
    using Field                     = VecField::field_type;
    constexpr static auto dimension = VecField::dimension;

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        if constexpr (dimension == 1)
            return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x);
        else if constexpr (dimension == 2)
            return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, vt_y, aL_y, aR_y, dL_y,
                                         dR_y);
        else if constexpr (dimension == 3)
            return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, vt_y, aL_y, aR_y, dL_y, dR_y,
                                         vt_z, aL_z, aR_z, dL_z, dR_z);
        else
            throw std::runtime_error(
                "Error - UpwindConstrainedTransportState - dimension not supported");
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        if constexpr (dimension == 1)
            return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x);
        else if constexpr (dimension == 2)
            return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, vt_y, aL_y, aR_y, dL_y,
                                         dR_y);
        else if constexpr (dimension == 3)
            return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, vt_y, aL_y, aR_y, dL_y, dR_y,
                                         vt_z, aL_z, aR_z, dL_z, dR_z);
        else
            throw std::runtime_error(
                "Error - UpwindConstrainedTransportState - dimension not supported");
    }

    template<auto direction>
    auto& getJt()
    {
        if constexpr (direction == Direction::X)
            return transverse().jt_x;
        else if constexpr (direction == Direction::Y)
            return transverse().jt_y;
        else if constexpr (direction == Direction::Z)
            return transverse().jt_z;
    }

    template<auto direction>
    auto& getRhot() const
    {
        if constexpr (direction == Direction::X)
            return transverse().rhot_x;
        else if constexpr (direction == Direction::Y)
            return transverse().rhot_y;
        else if constexpr (direction == Direction::Z)
            return transverse().rhot_z;
    }

    template<auto direction>
    void save(auto const vt, auto const& coefs, MeshIndex<dimension> const& idx)
    {
        auto assign_fields = [&](auto& vT, auto& aL, auto& aR, auto& dL, auto& dR) {
            vT(Component::X)(idx) = vt.x;
            vT(Component::Y)(idx) = vt.y;
            vT(Component::Z)(idx) = vt.z;

            aL(idx) = coefs[0];
            aR(idx) = coefs[1];
            dL(idx) = coefs[2];
            dR(idx) = coefs[3];
        };

        if constexpr (direction == Direction::X)
            assign_fields(vt_x, aL_x, aR_x, dL_x, dR_x);
        else if constexpr (direction == Direction::Y)
            assign_fields(vt_y, aL_y, aR_y, dL_y, dR_y);
        else if constexpr (direction == Direction::Z)
            assign_fields(vt_z, aL_z, aR_z, dL_z, dR_z);
    }

    template<auto direction>
    void save(auto const& vt, auto const& jt, auto const rhot, auto const& coefs,
              MeshIndex<dimension> const& idx)
    {
        save<direction>(vt, coefs, idx); // vt + uct coefficients

        auto& jT              = getJt<direction>();
        jT(Component::X)(idx) = jt.x;
        jT(Component::Y)(idx) = jt.y;
        jT(Component::Z)(idx) = jt.z;

        if constexpr (direction == Direction::X)
            transverse().rhot_x(idx) = rhot;
        else if constexpr (direction == Direction::Y)
            transverse().rhot_y(idx) = rhot;
        else if constexpr (direction == Direction::Z)
            transverse().rhot_z(idx) = rhot;
    }

    /// those 'transverse' velocity fields are an average value on faces (hence _x,
    /// _y, _z) of left and right states. The averaging formula depends on the Riemann solver.
    VecField vt_x{"v_t_x", MHDQuantity::Vector::VecFlux_x};
    VecField vt_y{"v_t_y", MHDQuantity::Vector::VecFlux_y};
    VecField vt_z{"v_t_z", MHDQuantity::Vector::VecFlux_z};

    /// advective and diffusive coefficients of the Riemann solver at the faces.
    /// a = advective / d = diffusive
    /// L = left / R = right
    /// _x = x face / _y = y face / _z = z face
    Field aL_x{"aL_x", MHDQuantity::Scalar::ScalarFlux_x},
        aR_x{"aR_x", MHDQuantity::Scalar::ScalarFlux_x},
        dL_x{"dL_x", MHDQuantity::Scalar::ScalarFlux_x},
        dR_x{"dR_x", MHDQuantity::Scalar::ScalarFlux_x};

    Field aL_y{"aL_y", MHDQuantity::Scalar::ScalarFlux_y},
        aR_y{"aR_y", MHDQuantity::Scalar::ScalarFlux_y},
        dL_y{"dL_y", MHDQuantity::Scalar::ScalarFlux_y},
        dR_y{"dR_y", MHDQuantity::Scalar::ScalarFlux_y};

    Field aL_z{"aL_z", MHDQuantity::Scalar::ScalarFlux_z},
        aR_z{"aR_z", MHDQuantity::Scalar::ScalarFlux_z},
        dL_z{"dL_z", MHDQuantity::Scalar::ScalarFlux_z},
        dR_z{"dR_z", MHDQuantity::Scalar::ScalarFlux_z};

private:
    auto& transverse() { return static_cast<Derived&>(*this).transverse_state(); }
    auto& transverse() const { return static_cast<Derived const&>(*this).transverse_state(); }
};

} // namespace PHARE::core


namespace PHARE::core::basic
{
// view of an UpwindConstrainedTransportTransverseState
template<typename VecField>
struct UpwindConstrainedTransportTransverseState
{
    using Field = VecField::field_type;

    VecField jt_x{}, jt_y{}, jt_z{};
    Field rhot_x{}, rhot_y{}, rhot_z{};
};

// view of an UpwindConstrainedTransportState (e.g. of one tile), the transverse state is null
// views if not Hall
template<typename VecField>
struct UpwindConstrainedTransportState
    : core::UpwindConstrainedTransportFields<UpwindConstrainedTransportState<VecField>, VecField>
{
    using transverse_type = UpwindConstrainedTransportTransverseState<VecField>;

    NO_DISCARD auto& transverse_state() { return transverse_; }
    NO_DISCARD auto& transverse_state() const { return transverse_; }

    transverse_type transverse_{};
    bool hall = false;
};

} // namespace PHARE::core::basic


namespace PHARE::core
{

template<typename VecField>
class UpwindConstrainedTransportState
    : public UpwindConstrainedTransportFields<UpwindConstrainedTransportState<VecField>, VecField>
{
    using Super
        = UpwindConstrainedTransportFields<UpwindConstrainedTransportState<VecField>, VecField>;
    using TransverseResistiveState = UpwindConstrainedTransportTransverseState<VecField>;

    friend Super;

public:
    using vecfield_type = VecField;

    template<typename VF>
    using view_t = basic::UpwindConstrainedTransportState<VF>;

    UpwindConstrainedTransportState() = default;
    explicit UpwindConstrainedTransportState(bool const isHall)
    {
        if (isHall)
            transverse_state_.emplace_back();
    }

    NO_DISCARD std::vector<TransverseResistiveState>& getRunTimeResourcesViewList()
    {
        return transverse_state_;
    }
    NO_DISCARD std::vector<TransverseResistiveState> const& getRunTimeResourcesViewList() const
    {
        return transverse_state_;
    }

    template<typename View>
    auto as(auto&& a, auto&&... args) const
    {
        using Transverse_v = View::transverse_type;

        auto const& s = static_cast<Super const&>(*this);
        return View{{{a(s.vt_x, args...)},
                     {a(s.vt_y, args...)},
                     {a(s.vt_z, args...)},
                     {a(s.aL_x, args...)},
                     {a(s.aR_x, args...)},
                     {a(s.dL_x, args...)},
                     {a(s.dR_x, args...)},
                     {a(s.aL_y, args...)},
                     {a(s.aR_y, args...)},
                     {a(s.dL_y, args...)},
                     {a(s.dR_y, args...)},
                     {a(s.aL_z, args...)},
                     {a(s.aR_z, args...)},
                     {a(s.dL_z, args...)},
                     {a(s.dR_z, args...)}},
                    transverse_state_.empty()
                        ? Transverse_v{}
                        : transverse_state_[0].template as<Transverse_v>(a, args...),
                    !transverse_state_.empty()};
    }

private:
    NO_DISCARD auto& transverse_state() { return transverse_state_[0]; }
    NO_DISCARD auto& transverse_state() const { return transverse_state_[0]; }

    // this vector is expected to hold at most one element. This is weird, but this is the only way
    // the resource manager can handle heterogeneous runtime resources for now.
    std::vector<TransverseResistiveState> transverse_state_;
};

} // namespace PHARE::core

#endif
