
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/resources_manager/resources_manager.hpp"

namespace PHARE
{
namespace amr
{
    /**
     * @brief offsetIsZero_ returns true of the transformation has zero offset
     */
    bool offsetIsZero(SAMRAI::hier::Transformation const& transformation)
    {
        auto const& offset = transformation.getOffset();
        auto dimension     = offset.getDim();
        return transformation.getOffset() == SAMRAI::hier::IntVector::getZero(dimension);
    }



    /**
     * @brief isSameBlock returns true if the transformation is not changing block
     */
    bool isSameBlock(SAMRAI::hier::Transformation const& transformation)
    {
        return transformation.getBeginBlock() == transformation.getEndBlock();
    }




    /**
     * @brief AMRToLocal sets the AMRBox to local indexing relative to the referenceAMRBox
     */
    SAMRAI::hier::Box& AMRToLocal(SAMRAI::hier::Box& AMRBox,
                                  SAMRAI::hier::Box const& referenceAMRBox)
    {
        AMRBox.setLower(AMRBox.lower() - referenceAMRBox.lower());
        AMRBox.setUpper(AMRBox.upper() - referenceAMRBox.lower());
        return AMRBox;
    }



    /**
     * @brief AMRToLocal returns a local indexed box relative to the referenceAMRBox from the AMRBox
     */
    SAMRAI::hier::Box AMRToLocal(SAMRAI::hier::Box const& AMRBox,
                                 SAMRAI::hier::Box const& referenceAMRBox)
    {
        SAMRAI::hier::Box localBox{AMRBox};
        localBox.setLower(AMRBox.lower() - referenceAMRBox.lower());
        localBox.setUpper(AMRBox.upper() - referenceAMRBox.lower());
        return localBox;
    }



    /**
     * @brief AMRToLocal returns the vector to add to a box to put it in the local index space
     * relative to the referenceAMRBox
     */
    SAMRAI::hier::IntVector AMRToLocal(SAMRAI::hier::Box const& referenceAMRBox)
    {
        SAMRAI::hier::Index zero{referenceAMRBox.getDim(), 0};
        return SAMRAI::hier::IntVector{zero - referenceAMRBox.lower()};
    }




    /**
     * @brief localToAMR returns the vector to add to a box to put it in AMR index space from a
     * local index relative to referenceAMRBox
     */
    SAMRAI::hier::IntVector localToAMRVector(SAMRAI::hier::Box const& referenceAMRBox)
    {
        return SAMRAI::hier::IntVector{referenceAMRBox.lower()};
    }
} // namespace amr


} // namespace PHARE


namespace PHARE::amr
{



ResourcesManagerGlobals& ResourcesManagerGlobals::INSTANCE()
{
    static ResourcesManagerGlobals globals;
    return globals;
}
} // namespace PHARE::amr
