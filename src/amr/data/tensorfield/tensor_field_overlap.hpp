#ifndef PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_OVERLAP_HPP
#define PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_OVERLAP_HPP


#include "core/data/tensorfield/tensorfield.hpp"
#include "amr/data/field/field_overlap.hpp"
#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/Transformation.h>

namespace PHARE::amr
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
template<std::size_t rank_ = 1>
class TensorFieldOverlap : public SAMRAI::hier::BoxOverlap
{
protected:
    auto constexpr static N = core::detail::tensor_field_dim_from_rank<rank_>();

public:
    static constexpr std::size_t rank = rank_;

    TensorFieldOverlap(std::array<std::shared_ptr<FieldOverlap>, N>&& overlaps)
        : transformation_{overlaps[0]->getTransformation()}
        , isOverlapEmpty_{true}
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            auto const& t = overlaps[i]->getTransformation();
            if (!transformations_equal_(t, transformation_))
            {
                throw std::runtime_error(
                    "Inconsistent transformation across FieldOverlap components.");
            }

            components_[i] = std::move(overlaps[i]);
            isOverlapEmpty_ &= components_[i]->isOverlapEmpty();
        }
    }

    ~TensorFieldOverlap() = default;



    bool isOverlapEmpty() const final { return isOverlapEmpty_; }



    const SAMRAI::hier::IntVector& getSourceOffset() const final
    {
        return transformation_.getOffset();
    }



    const SAMRAI::hier::Transformation& getTransformation() const final { return transformation_; }

    NO_DISCARD auto& operator[](std::size_t i) { return components_[i]; }
    NO_DISCARD auto& operator[](std::size_t i) const { return components_[i]; }

private:
    auto static _get_index_for(core::Component component)
    {
        auto val = static_cast<std::underlying_type_t<core::Component>>(component);
        if constexpr (rank == 1)
            return val;
        else if constexpr (rank == 2)
            return val - core::detail::tensor_field_dim_from_rank<1>();
    }

    bool transformations_equal_(const SAMRAI::hier::Transformation& a,
                                const SAMRAI::hier::Transformation& b)
    {
        return a.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE
               && b.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE
               && a.getOffset() == b.getOffset() && a.getBeginBlock() == b.getBeginBlock()
               && a.getEndBlock() == b.getEndBlock();
    }

    SAMRAI::hier::Transformation const transformation_;
    bool isOverlapEmpty_;

    std::array<std::shared_ptr<FieldOverlap>, N> components_;
};



} // namespace PHARE::amr

#endif
