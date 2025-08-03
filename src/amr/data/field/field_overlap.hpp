#ifndef PHARE_SRC_AMR_FIELD_FIELD_OVERLAP_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_OVERLAP_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/Transformation.h>

namespace PHARE
{
namespace amr
{
    /** \brief FieldOverlap is used to represent a region where data will be communicated betwen two
     * AMR patches
     *
     *  It will contain the exact form of the overlap between two patch for a fieldData with the
     * same quantity. It will also store any transformation between a source and destination patch.
     */
    /**
     * @brief The FieldOverlap class
     */
    class FieldOverlap : public SAMRAI::hier::BoxOverlap
    {
    public:
        FieldOverlap(SAMRAI::hier::BoxContainer const& boxes,
                     SAMRAI::hier::Transformation const& transformation)
            : destinationBoxes_{boxes}
            , transformation_{transformation}
            , isOverlapEmpty_{boxes.empty()}
        {
        }

        ~FieldOverlap() = default;



        bool isOverlapEmpty() const final { return isOverlapEmpty_; }



        const SAMRAI::hier::IntVector& getSourceOffset() const final
        {
            return transformation_.getOffset();
        }



        const SAMRAI::hier::Transformation& getTransformation() const final
        {
            return transformation_;
        }



        const SAMRAI::hier::BoxContainer& getDestinationBoxContainer() const
        {
            return destinationBoxes_;
        }


    private:
        SAMRAI::hier::BoxContainer const destinationBoxes_;
        SAMRAI::hier::Transformation const transformation_;
        bool const isOverlapEmpty_;
    };

} // namespace amr


} // namespace PHARE

#endif
