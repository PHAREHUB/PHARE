#ifndef PHARE_UPWIND_CONSTRAINED_TRANSPORT_UTILS_HPP
#define PHARE_UPWIND_CONSTRAINED_TRANSPORT_UTILS_HPP

#include "core/def.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include <stdexcept>

namespace PHARE::core
{

template<typename VecField>
class UpwindConstrainedTransportState
{
    using Field                     = VecField::field_type;
    constexpr static auto dimension = VecField::dimension;

public:
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
};

} // namespace PHARE::core

#endif
