#ifndef PHARE_UPWIND_CONSTRAINED_TRANSPORT_UTILS_HPP
#define PHARE_UPWIND_CONSTRAINED_TRANSPORT_UTILS_HPP

#include "core/def.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include <vector>

namespace PHARE::core
{

template<typename VecField>
class UpwindConstrainedTransportState
{
    using Field                     = VecField::field_type;
    constexpr static auto dimension = VecField::dimension;

    // jt (transverse interface current) and rhot (transverse interface density) are only needed
    // when the Hall term or resistivity is active. They are bundled in a holder so they can be
    // registered / allocated only when necessary, through a runtime resource list: the list is a
    // reference to persistent storage (this vector), so the resource manager sets buffers on these
    // members and not on temporaries. The manager recurses into the holder's compile-time list, so
    // its heterogeneous VecField + Field content is handled the same way as an ion population.
    struct TransverseResistiveState
    {
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

        VecField jt_x{"j_t_x", MHDQuantity::Vector::VecFlux_x};
        VecField jt_y{"j_t_y", MHDQuantity::Vector::VecFlux_y};
        VecField jt_z{"j_t_z", MHDQuantity::Vector::VecFlux_z};

        Field rhot_x{"rho_t_x", MHDQuantity::Scalar::ScalarFlux_x};
        Field rhot_y{"rho_t_y", MHDQuantity::Scalar::ScalarFlux_y};
        Field rhot_z{"rho_t_z", MHDQuantity::Scalar::ScalarFlux_z};
    };

public:
    UpwindConstrainedTransportState() = default;
    UpwindConstrainedTransportState(bool const isHall, bool const isResistive)
    {
        // hyper-resistivity implies the Hall term, so isHall || isResistive covers every case in
        // which jt / rhot are consumed (Hall EMF, resistive and hyper-resistive energy fluxes).
        if (isHall || isResistive)
            resistive_.emplace_back();
    }

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

    NO_DISCARD std::vector<TransverseResistiveState>& getRunTimeResourcesViewList()
    {
        return resistive_;
    }
    NO_DISCARD std::vector<TransverseResistiveState> const& getRunTimeResourcesViewList() const
    {
        return resistive_;
    }

    template<auto direction>
    auto& getJt()
    {
        if constexpr (direction == Direction::X)
            return resistive_[0].jt_x;
        else if constexpr (direction == Direction::Y)
            return resistive_[0].jt_y;
        else if constexpr (direction == Direction::Z)
            return resistive_[0].jt_z;
    }

    template<auto direction>
    auto& getRhot() const
    {
        if constexpr (direction == Direction::X)
            return resistive_[0].rhot_x;
        else if constexpr (direction == Direction::Y)
            return resistive_[0].rhot_y;
        else if constexpr (direction == Direction::Z)
            return resistive_[0].rhot_z;
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
            resistive_[0].rhot_x(idx) = rhot;
        else if constexpr (direction == Direction::Y)
            resistive_[0].rhot_y(idx) = rhot;
        else if constexpr (direction == Direction::Z)
            resistive_[0].rhot_z(idx) = rhot;
    }

    VecField vt_x{"v_t_x", MHDQuantity::Vector::VecFlux_x};
    VecField vt_y{"v_t_y", MHDQuantity::Vector::VecFlux_y};
    VecField vt_z{"v_t_z", MHDQuantity::Vector::VecFlux_z};

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
    std::vector<TransverseResistiveState> resistive_;
};

} // namespace PHARE::core

#endif
