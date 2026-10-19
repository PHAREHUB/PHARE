#ifndef PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP
#define PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP

#include "core/def.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/enum_table.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/mhd_quantities.hpp"


#include <array>

namespace PHARE::core::mhd
{

/**
 * @brief GridLayoutNdArrayImplYee class is a concrete GridLayoutNdArrayImpl used a Yee
 * type grid layout is needed.
 *
 * It provides methods related to grid layout operations:
 * - physical domain start/end indexes
 * - indexes of the first and last ghost nodes
 * - allocation sizes for Field attributes of other classes
 * - partial derivative operator (Faraday)
 * - physical coordinate given a field and a primal point (ix, iy, iz)
 * - cell centered coordinate given a primal point (ix, iy, iz)
 */
template<auto options>
class GridLayoutImplYee
{
    using Options = decltype(options);
    using Scalar  = Options::Scalar;

    static constexpr std::size_t dim = options.dimension;

public:
    struct GridData;
    static constexpr std::size_t dimension = dim;
    static constexpr std::string_view type = "yee";
    // The MHD layout reserves ghosts based on the reconstruction stencil width, plus the one
    // layer ampere loses by computing J on the shrinked ghost box.
    // static constexpr std::uint32_t reconstruction_nghosts = reconstruction_nghosts_;
    // Ghost width computed directly based on reconstruction stencil
    // static constexpr std::uint32_t ghost_width
    //     = nbrGhostsFromReconstruction<field_options.reconstruction_nghosts>();

    /**
     * @brief GridLayoutImpl<Selector<Layout,Layout::Yee>,dim>::initLayoutCentering_ initialize
     * the table MHDQuantityCentering_. This is THE important array in the GridLayout module.
     * This table knows which quantity is primal/dual along each direction. It is **this** array
     * that
     * **defines** what a Yee Layout is. Once this array is defined, the rest of the GridLayout
     * needs this array OK and can go on from here... hence all other functions in the Yee
     * interface are just calling private implementation common to all layouts
     */
    constexpr auto static initLayoutCentering_()
    {
        GridData const data{};

        std::array<QtyCentering, NBR_COMPO> const Rho = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Vx = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Vy = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Vz = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Bx = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const By = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Bz = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const B1x = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const B1y = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const B1z = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const P = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const rhoVx = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const rhoVy = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const rhoVz = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Etot  = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Etot1 = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Ex = {{data.dual, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Ey = {{data.primal, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Ez = {{data.primal, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Jx = {{data.dual, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Jy = {{data.primal, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Jz = {{data.primal, data.primal, data.dual}};


        std::array<QtyCentering, NBR_COMPO> const ScalarFlux_x
            = {{data.primal, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const ScalarFlux_y
            = {{data.dual, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const ScalarFlux_z
            = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const VecFluxX_x
            = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxY_x
            = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxZ_x
            = {{data.primal, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const VecFluxX_y
            = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxY_y
            = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxZ_y
            = {{data.dual, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const VecFluxX_z
            = {{data.dual, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxY_z
            = {{data.dual, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxZ_z
            = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const divB = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const ScalarAllPrimal
            = {{data.primal, data.primal, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const VecAllPrimalX
            = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const VecAllPrimalY
            = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const VecAllPrimalZ
            = {{data.primal, data.primal, data.primal}};

        using Qty                = MHDQuantity::Scalar;
        auto const _QtyCentering = make_enum_table<Qty, std::array<QtyCentering, NBR_COMPO>>({
            {Qty::rho, Rho},
            {Qty::Vx, Vx},
            {Qty::Vy, Vy},
            {Qty::Vz, Vz},
            {Qty::Bx, Bx},
            {Qty::By, By},
            {Qty::Bz, Bz},
            {Qty::B1x, B1x},
            {Qty::B1y, B1y},
            {Qty::B1z, B1z},
            {Qty::P, P},
            {Qty::Etot, Etot},
            {Qty::Etot1, Etot1},
            {Qty::rhoVx, rhoVx},
            {Qty::rhoVy, rhoVy},
            {Qty::rhoVz, rhoVz},
            {Qty::Ex, Ex},
            {Qty::Ey, Ey},
            {Qty::Ez, Ez},
            {Qty::Jx, Jx},
            {Qty::Jy, Jy},
            {Qty::Jz, Jz},
            {Qty::ScalarFlux_x, ScalarFlux_x},
            {Qty::ScalarFlux_y, ScalarFlux_y},
            {Qty::ScalarFlux_z, ScalarFlux_z},
            {Qty::VecFluxX_x, VecFluxX_x},
            {Qty::VecFluxY_x, VecFluxY_x},
            {Qty::VecFluxZ_x, VecFluxZ_x},
            {Qty::VecFluxX_y, VecFluxX_y},
            {Qty::VecFluxY_y, VecFluxY_y},
            {Qty::VecFluxZ_y, VecFluxZ_y},
            {Qty::VecFluxX_z, VecFluxX_z},
            {Qty::VecFluxY_z, VecFluxY_z},
            {Qty::VecFluxZ_z, VecFluxZ_z},
            {Qty::divB, divB},
            {Qty::ScalarAllPrimal, ScalarAllPrimal},
            {Qty::VecAllPrimalX, VecAllPrimalX},
            {Qty::VecAllPrimalY, VecAllPrimalY},
            {Qty::VecAllPrimalZ, VecAllPrimalZ},
        });

        return _QtyCentering;
    }

    //! says for each MHDQuantity::Quantity whether it is primal or dual, in each direction
    constexpr static std::array<std::array<QtyCentering, NBR_COMPO>,
                                static_cast<std::size_t>(MHDQuantity::Scalar::count)> const
        qtyCentering_{initLayoutCentering_()};

    static std::size_t const dim_{dim};

    //! qtyCentering_ restricted to the first dim directions, indexed by MHDQuantity::Scalar
    constexpr static auto scalarCentering_ = [] {
        std::array<std::array<QtyCentering, dim>, qtyCentering_.size()> centerings{};
        for (std::size_t iQty = 0; iQty < centerings.size(); ++iQty)
            for (std::size_t iDir = 0; iDir < dim; ++iDir)
                centerings[iQty][iDir] = qtyCentering_[iQty][iDir];
        return centerings;
    }();

    //! centerings of each vector's components, indexed by MHDQuantity::Vector
    constexpr static auto vectorCentering_ = [] {
        std::array<std::array<std::array<QtyCentering, dim>, 3>,
                   MHDQuantity::vector_components.size()>
            centerings{};
        for (std::size_t iVec = 0; iVec < centerings.size(); ++iVec)
            for (std::size_t iComp = 0; iComp < 3; ++iComp)
                centerings[iVec][iComp] = scalarCentering_[static_cast<std::size_t>(
                    MHDQuantity::vector_components[iVec][iComp])];
        return centerings;
    }();

    // ------------------------------------------------------------------------
    //                          PUBLIC INTERFACE
    // ------------------------------------------------------------------------
public:
    NO_DISCARD constexpr static std::array<QtyCentering, dim>
    centering(MHDQuantity::Scalar const qty)
    {
        auto const idx = static_cast<std::size_t>(qty);
        if (idx >= scalarCentering_.size())
            throw std::runtime_error("Wrong MHDQuantity");
        return scalarCentering_[idx];
    }

    NO_DISCARD constexpr static std::array<std::array<QtyCentering, dim>, 3>
    centering(MHDQuantity::Vector const qty)
    {
        auto const idx = static_cast<std::size_t>(qty);
        if (idx >= vectorCentering_.size())
            throw std::runtime_error("Wrong MHDQuantity");
        return vectorCentering_[idx];
    }
};



template<auto options>
struct GridLayoutImplYee<options>::GridData
{
    static constexpr Direction dirX = Direction::X;
    static constexpr Direction dirY = Direction::Y;
    static constexpr Direction dirZ = Direction::Z;

    static constexpr QtyCentering primal = QtyCentering::primal;
    static constexpr QtyCentering dual   = QtyCentering::dual;

    static constexpr std::uint32_t idirX = static_cast<std::uint32_t>(Direction::X);
    static constexpr std::uint32_t idirY = static_cast<std::uint32_t>(Direction::Y);
    static constexpr std::uint32_t idirZ = static_cast<std::uint32_t>(Direction::Z);

    static constexpr std::uint32_t irho = static_cast<std::uint32_t>(MHDQuantity::Scalar::rho);

    static constexpr std::uint32_t iVx = static_cast<std::uint32_t>(MHDQuantity::Scalar::Vx);
    static constexpr std::uint32_t iVy = static_cast<std::uint32_t>(MHDQuantity::Scalar::Vy);
    static constexpr std::uint32_t iVz = static_cast<std::uint32_t>(MHDQuantity::Scalar::Vz);

    static constexpr std::uint32_t iBx = static_cast<std::uint32_t>(MHDQuantity::Scalar::Bx);
    static constexpr std::uint32_t iBy = static_cast<std::uint32_t>(MHDQuantity::Scalar::By);
    static constexpr std::uint32_t iBz = static_cast<std::uint32_t>(MHDQuantity::Scalar::Bz);

    static constexpr std::uint32_t iB1x = static_cast<std::uint32_t>(MHDQuantity::Scalar::B1x);
    static constexpr std::uint32_t iB1y = static_cast<std::uint32_t>(MHDQuantity::Scalar::B1y);
    static constexpr std::uint32_t iB1z = static_cast<std::uint32_t>(MHDQuantity::Scalar::B1z);

    static constexpr std::uint32_t iP = static_cast<std::uint32_t>(MHDQuantity::Scalar::P);

    static constexpr std::uint32_t iEtot = static_cast<std::uint32_t>(MHDQuantity::Scalar::Etot);
    static constexpr std::uint32_t iEtot1 = static_cast<std::uint32_t>(MHDQuantity::Scalar::Etot1);

    static constexpr std::uint32_t irhoVx = static_cast<std::uint32_t>(MHDQuantity::Scalar::rhoVx);
    static constexpr std::uint32_t irhoVy = static_cast<std::uint32_t>(MHDQuantity::Scalar::rhoVy);
    static constexpr std::uint32_t irhoVz = static_cast<std::uint32_t>(MHDQuantity::Scalar::rhoVz);

    static constexpr std::uint32_t iEx = static_cast<std::uint32_t>(MHDQuantity::Scalar::Ex);
    static constexpr std::uint32_t iEy = static_cast<std::uint32_t>(MHDQuantity::Scalar::Ey);
    static constexpr std::uint32_t iEz = static_cast<std::uint32_t>(MHDQuantity::Scalar::Ez);

    static constexpr std::uint32_t iJx = static_cast<std::uint32_t>(MHDQuantity::Scalar::Jx);
    static constexpr std::uint32_t iJy = static_cast<std::uint32_t>(MHDQuantity::Scalar::Jy);
    static constexpr std::uint32_t iJz = static_cast<std::uint32_t>(MHDQuantity::Scalar::Jz);

    static constexpr std::uint32_t iScalarFlux_x
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::ScalarFlux_x);
    static constexpr std::uint32_t iScalarFlux_y
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::ScalarFlux_y);
    static constexpr std::uint32_t iScalarFlux_z
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::ScalarFlux_z);

    static constexpr std::uint32_t iVecFluxX_x
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxX_x);
    static constexpr std::uint32_t iVecFluxY_x
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxY_x);
    static constexpr std::uint32_t iVecFluxZ_x
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxZ_x);

    static constexpr std::uint32_t iVecFluxX_y
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxX_y);
    static constexpr std::uint32_t iVecFluxY_y
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxY_y);
    static constexpr std::uint32_t iVecFluxZ_y
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxZ_y);

    static constexpr std::uint32_t iVecFluxX_z
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxX_z);
    static constexpr std::uint32_t iVecFluxY_z
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxY_z);
    static constexpr std::uint32_t iVecFluxZ_z
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxZ_z);

    static constexpr std::uint32_t idivB = static_cast<std::uint32_t>(MHDQuantity::Scalar::divB);

    static constexpr std::uint32_t iScalarAllPrimal
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::ScalarAllPrimal);

    static constexpr std::uint32_t iVecAllPrimalX
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecAllPrimalX);
    static constexpr std::uint32_t iVecAllPrimalY
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecAllPrimalY);
    static constexpr std::uint32_t iVecAllPrimalZ
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecAllPrimalZ);
};



} // namespace PHARE::core::mhd


#endif // PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP
