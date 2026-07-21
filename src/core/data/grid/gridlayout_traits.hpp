#ifndef PHARE_CORE_GRID_GRIDLAYOUT_TRAITS_HPP
#define PHARE_CORE_GRID_GRIDLAYOUT_TRAITS_HPP

#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include <array>
#include <tuple>
#include <concepts>

namespace PHARE::core
{
/**
 * @brief Define the requirements for a grid layout type.
 *
 * @see GridLayout
 *
 */
template<typename GridLayoutT>
concept IsGridLayout
    = requires(GridLayoutT g, Direction d, QtyCentering c, int idx, std::uint32_t u_idx,
               typename GridLayoutT::Quantity::Scalar qty,
               typename GridLayoutT::Quantity::Vector v_qty,
               Point<int, GridLayoutT::dimension> p_int, Box<int, GridLayoutT::dimension> b_int) {
          { GridLayoutT::dimension } -> std::convertible_to<std::size_t>;
          { GridLayoutT::interp_order } -> std::convertible_to<std::size_t>;
          { GridLayoutT::nbrParticleGhosts() } -> std::convertible_to<std::size_t>;

          { g.origin() } -> std::same_as<Point<double, GridLayoutT::dimension>>;
          { g.meshSize() } -> std::same_as<std::array<double, GridLayoutT::dimension> const&>;
          { g.inverseMeshSize(d) } -> std::convertible_to<double>;
          { g.inverseMeshSize() } -> std::same_as<std::array<double, GridLayoutT::dimension>>;
          {
              g.nbrCells()
          } -> std::convertible_to<std::array<std::uint32_t, GridLayoutT::dimension> const&>;
          { g.cellVolume() } -> std::convertible_to<double>;
          { g.layoutName() } -> std::same_as<std::string>;
          { g.levelNumber() } -> std::convertible_to<int>;

          { g.AMRBox() } -> std::same_as<Box<int, GridLayoutT::dimension> const&>;
          { g.localToAMR(p_int) } -> std::same_as<Point<int, GridLayoutT::dimension>>;
          { g.localToAMR(b_int) } -> std::same_as<Box<int, GridLayoutT::dimension>>;
          { g.AMRToLocal(p_int) } -> std::same_as<Point<std::uint32_t, GridLayoutT::dimension>>;
          { g.AMRToLocal(b_int) } -> std::same_as<Box<std::uint32_t, GridLayoutT::dimension>>;

          {
              GridLayoutT::centering(qty)
          } -> std::same_as<std::array<QtyCentering, GridLayoutT::dimension>>;
          {
              GridLayoutT::centering(v_qty)
          } -> std::same_as<std::array<std::array<QtyCentering, GridLayoutT::dimension>, 3>>;
          { GridLayoutT::changeCentering(c) } -> std::same_as<QtyCentering>;
          { GridLayoutT::nextIndex(c, u_idx) } -> std::convertible_to<std::uint32_t>;
          { GridLayoutT::prevIndex(c, u_idx) } -> std::convertible_to<std::uint32_t>;

          { g.physicalStartIndex(c, d) } -> std::convertible_to<std::uint32_t>;
          { g.physicalEndIndex(c, d) } -> std::convertible_to<std::uint32_t>;
          { g.ghostStartIndex(c, d) } -> std::convertible_to<std::uint32_t>;
          { g.ghostEndIndex(c, d) } -> std::convertible_to<std::uint32_t>;
          { g.ghostStartToEnd(c, d) } -> std::same_as<std::tuple<std::uint32_t, std::uint32_t>>;
          { g.physicalStartToEnd(c, d) } -> std::same_as<std::tuple<std::uint32_t, std::uint32_t>>;

          { g.allocSize(qty) } -> std::same_as<std::array<std::uint32_t, GridLayoutT::dimension>>;
          {
              g.allocSizeDerived(qty, d)
          } -> std::same_as<std::array<std::uint32_t, GridLayoutT::dimension>>;
          {
              g.nbrPhysicalNodes(qty)
          } -> std::same_as<std::array<std::uint32_t, GridLayoutT::dimension>>;
          { GridLayoutT::nbrGhosts(c) } -> std::convertible_to<std::uint32_t>;

          {
              g.cellCenteredCoordinates(p_int)
          } -> std::same_as<Point<double, GridLayoutT::dimension>>;

          // --- 8. Projection Helpers (Static Methods) ---
          //   { GridLayoutT::momentsToEx() };
          //   { GridLayoutT::momentsToEy() };
          //   { GridLayoutT::momentsToEz() };
          //   { GridLayoutT::ExToMoments() };
          //   { GridLayoutT::EyToMoments() };
          //   { GridLayoutT::EzToMoments() };
          //   { GridLayoutT::JxToMoments() };
          //   { GridLayoutT::JyToMoments() };
          //   { GridLayoutT::JzToMoments() };
          //   { GridLayoutT::BxToEx() };
          //   { GridLayoutT::ByToEx() };
          //   { GridLayoutT::BzToEx() };
          //   { GridLayoutT::JxToEx() };
          //   { GridLayoutT::JyToEy() };
          //   { GridLayoutT::JzToEz() };
          //   { GridLayoutT::faceXToCellCenter() };
          //   { GridLayoutT::edgeXToCellCenter() };

          { g.amr_lcl_idx(b_int) };
          { g.amr_lcl_idx() };
      };

} // namespace PHARE::core

#endif // PHARE_CORE_GRID_GRIDLAYOUT_TRAITS_HPP
