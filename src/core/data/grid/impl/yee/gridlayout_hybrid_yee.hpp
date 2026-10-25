#ifndef PHARE_CORE_GRID_IMPL_YEE_GRIDLAYOUT_HYBRID_YEE_HPP
#define PHARE_CORE_GRID_IMPL_YEE_GRIDLAYOUT_HYBRID_YEE_HPP


#include "core/def.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/enum_table.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"


#include <array>


namespace PHARE::core::hybrid
{


/**
 * @brief GridData provides constants used to initialize:
 * - hybridQuantity centerings
 * - physical start/end indexes
 * - ghost start/end indexes
 * - numbers of padding cells and physical cells
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


    /**
     * @brief GridLayoutImpl<Selector<Layout,Layout::Yee>,dim>::initLayoutCentering_ initialize
     * the table hybridQuantityCentering_. This is THE important array in the GridLayout module.
     * This table knows which quantity is primal/dual along each direction. It is **this** array
     * that
     * **defines** what a Yee Layout is. Once this array is defined, the rest of the GridLayout
     * needs this array OK and can go on from here... hence all other functions in the Yee
     * interface are just calling private implementation common to all layouts
     */
    constexpr auto static initLayoutCentering_()
    {
        GridData const data{};
        std::array<QtyCentering, NBR_COMPO> const Bx = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const By = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Bz = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const Ex = {{data.dual, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Ey = {{data.primal, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Ez = {{data.primal, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Jx = {{data.dual, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Jy = {{data.primal, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Jz = {{data.primal, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Rho = {{data.primal, data.primal, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const Vx = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Vy = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Vz = {{data.primal, data.primal, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const Mxx = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Mxy = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Mxz = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Myy = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Myz = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Mzz = {{data.primal, data.primal, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const P = {{data.primal, data.primal, data.primal}};

        using Qty                     = HybridQuantity::Scalar;
        auto const hybridQtyCentering = make_enum_table<Qty, std::array<QtyCentering, NBR_COMPO>>({
            {Qty::Bx, Bx},   {Qty::By, By},   {Qty::Bz, Bz},   {Qty::Ex, Ex},   {Qty::Ey, Ey},
            {Qty::Ez, Ez},   {Qty::Jx, Jx},   {Qty::Jy, Jy},   {Qty::Jz, Jz},   {Qty::rho, Rho},
            {Qty::Vx, Vx},   {Qty::Vy, Vy},   {Qty::Vz, Vz},   {Qty::P, P},     {Qty::Mxx, Mxx},
            {Qty::Mxy, Mxy}, {Qty::Mxz, Mxz}, {Qty::Myy, Myy}, {Qty::Myz, Myz}, {Qty::Mzz, Mzz},
        });

        return hybridQtyCentering;
    }



    //! says for each HybridQuantity::Quantity whether it is primal or dual, in each direction
    constexpr static std::array<std::array<QtyCentering, NBR_COMPO>,
                                static_cast<std::size_t>(HybridQuantity::Scalar::count)> const
        qtyCentering_{initLayoutCentering_()};

    static std::size_t const dim_{dim};

    //! qtyCentering_ restricted to the first dim directions, indexed by HybridQuantity::Scalar
    constexpr static auto scalarCentering_ = [] {
        std::array<std::array<QtyCentering, dim>, qtyCentering_.size()> centerings{};
        for (std::size_t iQty = 0; iQty < centerings.size(); ++iQty)
            for (std::size_t iDir = 0; iDir < dim; ++iDir)
                centerings[iQty][iDir] = qtyCentering_[iQty][iDir];
        return centerings;
    }();

    //! centerings of each vector's components, indexed by HybridQuantity::Vector
    constexpr static auto vectorCentering_ = [] {
        std::array<std::array<std::array<QtyCentering, dim>, 3>,
                   HybridQuantity::vector_components.size()>
            centerings{};
        for (std::size_t iVec = 0; iVec < centerings.size(); ++iVec)
            for (std::size_t iComp = 0; iComp < 3; ++iComp)
                centerings[iVec][iComp] = scalarCentering_[static_cast<std::size_t>(
                    HybridQuantity::vector_components[iVec][iComp])];
        return centerings;
    }();

    // ------------------------------------------------------------------------
    //                          PUBLIC INTERFACE
    // ------------------------------------------------------------------------
public:
    NO_DISCARD constexpr static std::array<QtyCentering, dim>
    centering(HybridQuantity::Scalar const qty)
    {
        auto const idx = static_cast<std::size_t>(qty);
        if (idx >= scalarCentering_.size())
            throw std::runtime_error("Wrong hybridQuantity");
        return scalarCentering_[idx];
    }

    NO_DISCARD constexpr static std::array<std::array<QtyCentering, dim>, 3>
    centering(HybridQuantity::Vector const qty)
    {
        auto const idx = static_cast<std::size_t>(qty);
        if (idx >= vectorCentering_.size())
            throw std::runtime_error("Wrong hybridQuantity");
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

    static constexpr std::uint32_t iBx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Bx);
    static constexpr std::uint32_t iBy = static_cast<std::uint32_t>(HybridQuantity::Scalar::By);
    static constexpr std::uint32_t iBz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Bz);

    static constexpr std::uint32_t iEx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ex);
    static constexpr std::uint32_t iEy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ey);
    static constexpr std::uint32_t iEz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ez);

    static constexpr std::uint32_t iJx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jx);
    static constexpr std::uint32_t iJy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jy);
    static constexpr std::uint32_t iJz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jz);

    static constexpr std::uint32_t irho = static_cast<std::uint32_t>(HybridQuantity::Scalar::rho);

    static constexpr std::uint32_t iVx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vx);
    static constexpr std::uint32_t iVy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vy);
    static constexpr std::uint32_t iVz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vz);

    static constexpr std::uint32_t iMxx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxx);
    static constexpr std::uint32_t iMxy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxy);
    static constexpr std::uint32_t iMxz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxz);
    static constexpr std::uint32_t iMyy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Myy);
    static constexpr std::uint32_t iMyz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Myz);
    static constexpr std::uint32_t iMzz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mzz);

    static constexpr std::uint32_t iP = static_cast<std::uint32_t>(HybridQuantity::Scalar::P);
};



} // namespace PHARE::core::hybrid


#endif // PHARE_CORE_GRID_IMPL_YEE_GRIDLAYOUT_HYBRID_YEE_HPP
