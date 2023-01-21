#ifndef PHARE_MAGNETIC_FIELD_COARSENER
#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto a
     * given coarse node
     *
     * A MagneticFieldCoarsener object is created each time the refine() method of the
     * FieldCoarsenOperator is called and its operator() is called for each coarse index.
     * It is the default coarsening policy and used for any field that does not come with
     * specific constraints (such as conserving some property in the coarsening process).
     */
    template<std::size_t dimension>
    class MagneticFieldCoarsener
    {
    public:
        MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                               SAMRAI::hier::Box const& sourceBox,
                               SAMRAI::hier::Box const& destinationBox,
                               SAMRAI::hier::IntVector const& ratio)
        {
        }
        template<typename FieldT>
        void operator()(FieldT const& fineField, FieldT& coarseField,
                        core::Point<int, dimension> coarseIndex)
        {
        }
    };
} // namespace amr
} // namespace PHARE
#endif
