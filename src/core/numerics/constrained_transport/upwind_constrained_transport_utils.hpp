#ifndef PHARE_UPWIND_CONSTRAINED_TRANSPORT_UTILS_HPP
#define PHARE_UPWIND_CONSTRAINED_TRANSPORT_UTILS_HPP

#include "core/def.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

namespace PHARE::core
{

template<typename VecField, bool Hall, bool Resistivity>
class UpwindConstrainedTransportState
{
    using Field                     = VecField::field_type;
    constexpr static auto dimension = VecField::dimension;

public:
    UpwindConstrainedTransportState() = default;

    // Edge-B buffers (Bt_*_at_E*) carry upwind-averaged transverse B at E-field edge
    // locations for the CT Poynting energy correction. Bt_z_at_Ey / Bt_y_at_Ez are needed
    // in all dimensions; Bt_x_at_Ez / Bt_z_at_Ex in 2D+; Bt_x_at_Ey / Bt_y_at_Ex in 3D.
    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        if constexpr (dimension == 1)
        {
            if constexpr (Hall || Resistivity)
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, jt_x, rhot_x, Bt_z_at_Ey,
                                             Bt_y_at_Ez);
            else
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, Bt_z_at_Ey, Bt_y_at_Ez);
        }
        else if constexpr (dimension == 2)
        {
            if constexpr (Hall || Resistivity)
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, jt_x, rhot_x, vt_y, aL_y,
                                             aR_y, dL_y, dR_y, jt_y, rhot_y, Bt_z_at_Ey, Bt_y_at_Ez,
                                             Bt_x_at_Ez, Bt_z_at_Ex);
            else
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, vt_y, aL_y, aR_y, dL_y,
                                             dR_y, Bt_z_at_Ey, Bt_y_at_Ez, Bt_x_at_Ez, Bt_z_at_Ex);
        }
        else if constexpr (dimension == 3)
        {
            if constexpr (Hall || Resistivity)
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, jt_x, rhot_x, vt_y, aL_y,
                                             aR_y, dL_y, dR_y, jt_y, rhot_y, vt_z, aL_z, aR_z, dL_z,
                                             dR_z, jt_z, rhot_z, Bt_z_at_Ey, Bt_y_at_Ez, Bt_x_at_Ez,
                                             Bt_z_at_Ex, Bt_x_at_Ey, Bt_y_at_Ex);
            else
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, vt_y, aL_y, aR_y, dL_y,
                                             dR_y, vt_z, aL_z, aR_z, dL_z, dR_z, Bt_z_at_Ey,
                                             Bt_y_at_Ez, Bt_x_at_Ez, Bt_z_at_Ex, Bt_x_at_Ey,
                                             Bt_y_at_Ex);
        }
        else
            throw std::runtime_error(
                "Error - UpwindConstrainedTransportState - dimension not supported");
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        if constexpr (dimension == 1)
        {
            if constexpr (Hall || Resistivity)
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, jt_x, rhot_x, Bt_z_at_Ey,
                                             Bt_y_at_Ez);
            else
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, Bt_z_at_Ey, Bt_y_at_Ez);
        }
        else if constexpr (dimension == 2)
        {
            if constexpr (Hall || Resistivity)
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, jt_x, rhot_x, vt_y, aL_y,
                                             aR_y, dL_y, dR_y, jt_y, rhot_y, Bt_z_at_Ey, Bt_y_at_Ez,
                                             Bt_x_at_Ez, Bt_z_at_Ex);
            else
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, vt_y, aL_y, aR_y, dL_y,
                                             dR_y, Bt_z_at_Ey, Bt_y_at_Ez, Bt_x_at_Ez, Bt_z_at_Ex);
        }
        else if constexpr (dimension == 3)
        {
            if constexpr (Hall || Resistivity)
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, jt_x, rhot_x, vt_y, aL_y,
                                             aR_y, dL_y, dR_y, jt_y, rhot_y, vt_z, aL_z, aR_z, dL_z,
                                             dR_z, jt_z, rhot_z, Bt_z_at_Ey, Bt_y_at_Ez, Bt_x_at_Ez,
                                             Bt_z_at_Ex, Bt_x_at_Ey, Bt_y_at_Ex);
            else
                return std::forward_as_tuple(vt_x, aL_x, aR_x, dL_x, dR_x, vt_y, aL_y, aR_y, dL_y,
                                             dR_y, vt_z, aL_z, aR_z, dL_z, dR_z, Bt_z_at_Ey,
                                             Bt_y_at_Ez, Bt_x_at_Ez, Bt_z_at_Ex, Bt_x_at_Ey,
                                             Bt_y_at_Ex);
        }
        else
            throw std::runtime_error(
                "Error - UpwindConstrainedTransportState - dimension not supported");
    }

    // Getters for edge-centered B fields (CT Poynting energy correction).
    auto const& getBx_at_Ez() const { return Bt_x_at_Ez; }
    auto const& getBy_at_Ez() const { return Bt_y_at_Ez; }
    auto const& getBz_at_Ey() const { return Bt_z_at_Ey; }
    auto const& getBx_at_Ey() const { return Bt_x_at_Ey; }
    auto const& getBz_at_Ex() const { return Bt_z_at_Ex; }
    auto const& getBy_at_Ex() const { return Bt_y_at_Ex; }

    template<auto direction>
    auto& getJt()
    {
        if constexpr (direction == Direction::X)
            return jt_x;
        else if constexpr (direction == Direction::Y)
            return jt_y;
        else if constexpr (direction == Direction::Z)
            return jt_z;
    }

    template<auto direction>
    auto& getRhot() const
    {
        if constexpr (direction == Direction::X)
            return rhot_x;
        else if constexpr (direction == Direction::Y)
            return rhot_y;
        else if constexpr (direction == Direction::Z)
            return rhot_z;
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
        auto assign_fields
            = [&](auto& vT, auto& jT, auto& rhoT, auto& aL, auto& aR, auto& dL, auto& dR) {
                  vT(Component::X)(idx) = vt.x;
                  vT(Component::Y)(idx) = vt.y;
                  vT(Component::Z)(idx) = vt.z;

                  jT(Component::X)(idx) = jt.x;
                  jT(Component::Y)(idx) = jt.y;
                  jT(Component::Z)(idx) = jt.z;

                  rhoT(idx) = rhot;

                  aL(idx) = coefs[0];
                  aR(idx) = coefs[1];
                  dL(idx) = coefs[2];
                  dR(idx) = coefs[3];
              };

        if constexpr (direction == Direction::X)
            assign_fields(vt_x, jt_x, rhot_x, aL_x, aR_x, dL_x, dR_x);
        else if constexpr (direction == Direction::Y)
            assign_fields(vt_y, jt_y, rhot_y, aL_y, aR_y, dL_y, dR_y);
        else if constexpr (direction == Direction::Z)
            assign_fields(vt_z, jt_z, rhot_z, aL_z, aR_z, dL_z, dR_z);
    }

    VecField vt_x{"v_t_x", MHDQuantity::Vector::VecFlux_x};
    VecField vt_y{"v_t_y", MHDQuantity::Vector::VecFlux_y};
    VecField vt_z{"v_t_z", MHDQuantity::Vector::VecFlux_z};

    VecField jt_x{"j_t_x", MHDQuantity::Vector::VecFlux_x};
    VecField jt_y{"j_t_y", MHDQuantity::Vector::VecFlux_y};
    VecField jt_z{"j_t_z", MHDQuantity::Vector::VecFlux_z};

    Field rhot_x{"rho_t_x", MHDQuantity::Scalar::ScalarFlux_x};
    Field rhot_y{"rho_t_y", MHDQuantity::Scalar::ScalarFlux_y};
    Field rhot_z{"rho_t_z", MHDQuantity::Scalar::ScalarFlux_z};

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

    // Edge-centered B fields for the Poynting energy correction, stored at the E-field
    // edge locations where they are computed during CT.
    Field Bt_x_at_Ez{"Bx1ez", MHDQuantity::Scalar::Ez}; // Bx at z-edge
    Field Bt_y_at_Ez{"Bx2ez", MHDQuantity::Scalar::Ez}; // By at z-edge
    Field Bt_z_at_Ey{"Bx3ey", MHDQuantity::Scalar::Ey}; // Bz at y-edge (x-face in 2D)
    Field Bt_z_at_Ex{"Bx3ex", MHDQuantity::Scalar::Ex}; // Bz at x-edge (y-face in 2D)
    Field Bt_x_at_Ey{"Bx1ey", MHDQuantity::Scalar::Ey}; // Bx at y-edge (3D only)
    Field Bt_y_at_Ex{"Bx2ex", MHDQuantity::Scalar::Ex}; // By at x-edge (3D only)
};

} // namespace PHARE::core

#endif
