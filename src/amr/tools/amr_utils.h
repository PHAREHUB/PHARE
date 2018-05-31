#ifndef UTILS_H
#define UTILS_H

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>



namespace PHARE
{
/**
 * @brief offsetIsZero_ returns true of the transformation has zero offset
 */
bool offsetIsZero(SAMRAI::hier::Transformation const& transformation);



/**
 * @brief isSameBlock returns true if the transformation is not changing block
 */
bool isSameBlock(SAMRAI::hier::Transformation const& transformation);



/**
 * @brief AMRToLocal sets the AMRBox to local indexing relative to the referenceAMRBox
 */
void AMRToLocal(SAMRAI::hier::Box& AMRBox, SAMRAI::hier::Box const& referenceAMRBox);



/**
 * @brief AMRToLocal returns a local indexed box relative to the referenceAMRBox from the AMRBox
 */
SAMRAI::hier::Box AMRToLocal(SAMRAI::hier::Box const& AMRBox,
                             SAMRAI::hier::Box const& referenceAMRBox);



/**
 * @brief AMRToLocal returns the vector to add to a box to put it in the local index space
 * relative to the referenceAMRBox
 */
SAMRAI::hier::IntVector AMRToLocal(SAMRAI::hier::Box const& referenceAMRBox);


/**
 * @brief localToAMR returns the vector to add to a box to put it in AMR index space from a
 * local index relative to referenceAMRBox
 */
SAMRAI::hier::IntVector localToAMR(SAMRAI::hier::Box const& referenceAMRBox);

} // namespace PHARE

#endif // UTILS_H
