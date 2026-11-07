#ifndef PHARE_AMR_MESSENGERS_FIELD_BORDER_SUM_HPP
#define PHARE_AMR_MESSENGERS_FIELD_BORDER_SUM_HPP

#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_overlap.hpp"
#include "amr/data/tensorfield/tensor_field_overlap.hpp"
#include "amr/messengers/field_operate_transaction.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/GlobalId.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/Transformation.h>
#include <SAMRAI/xfer/RefineClasses.h>
#include <SAMRAI/xfer/RefineTransactionFactory.h>

#include <map>
#include <span>
#include <memory>
#include <vector>
#include <cstddef>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <array>
#include <type_traits>

namespace PHARE::amr
{
/* Border sums whose result does not depend on the patch computing them.
 *
 * A node in the ghost boxes of several patches of a level holds a partial deposit on each of
 * them, and every patch holding it needs the sum of these deposits. Floating point addition is
 * not associative: added in the order the transactions of a schedule run, which differs from one
 * holder to another, the holders of a node with three contributions or more disagree in the last
 * bits. Here the contributions to a node are added in one order, the same on every holder: that
 * of the contributing patch boxes (lower, then upper corner, lexicographic), seen from the node.
 * Two holders see the same contributing boxes, up to one common periodic shift, which does not
 * change their order.
 *
 * The transactions of the schedule store the contribution they carry instead of adding it, and
 * BorderSums::fold, once the schedule has run, adds them up in that order on each patch, in place.
 * The order is set once per schedule, which is rebuilt at regrid.
 */


namespace detail
{
    // the components of an overlap of FieldData (one) or TensorFieldData (one per component)
    template<typename Overlap>
    auto overlapComponents(Overlap const& overlap)
    {
        if constexpr (std::is_same_v<Overlap, FieldOverlap>)
            return std::array<FieldOverlap const*, 1>{&overlap};
        else
            return core::for_N_make_array<
                core::detail::tensor_field_dim_from_rank<Overlap::rank>()>(
                [&](auto i) -> FieldOverlap const* { return overlap[i].get(); });
    }

    template<typename Overlap>
    auto makeOverlap(auto const& componentBoxes, SAMRAI::hier::Transformation const& transformation)
    {
        if constexpr (std::is_same_v<Overlap, FieldOverlap>)
            return std::make_shared<FieldOverlap>(componentBoxes[0], transformation);
        else
        {
            auto constexpr N = std::tuple_size_v<std::decay_t<decltype(componentBoxes)>>;
            std::array<std::shared_ptr<FieldOverlap>, N> components;
            for (std::size_t c = 0; c < N; ++c)
                components[c] = std::make_shared<FieldOverlap>(componentBoxes[c], transformation);
            return std::make_shared<Overlap>(std::move(components));
        }
    }

    // number of values a stream carries for `overlap`: one per node of each destination box
    std::size_t streamSize(auto const& overlap)
    {
        std::size_t size = 0;
        for (auto const* component : overlapComponents(overlap))
            size += component->getDestinationBoxContainer().getTotalSizeOfBoxes();
        return size;
    }

    // the overlap type of FieldData (FieldOverlap) or TensorFieldData
    template<typename FieldData_t>
    struct OverlapOf
    {
        using type = FieldOverlap;
    };
    template<typename FieldData_t>
        requires requires { FieldData_t::rank; }
    struct OverlapOf<FieldData_t>
    {
        using type = TensorFieldOverlap<FieldData_t::rank>;
    };

    template<typename T>
    struct SetZero
    {
        void operator()(auto const& /*ignored*/) { t = T{0}; }
        T& t;
    };
} // namespace detail



template<typename FieldData_t>
class BorderSums
{
public:
    using value_type = FieldData_t::value_type;
    using Overlap_t  = detail::OverlapOf<FieldData_t>::type;

private:
    auto constexpr static dimension = FieldData_t::dimension;
    using Box_t                     = core::Box<int, dimension>;

    struct Contribution
    {
        Box_t source; // contributing patch box, in the destination index space
        std::shared_ptr<Overlap_t const> overlap;
        std::size_t offset; // of its values in the store
        std::size_t size;
    };

    struct Destination
    {
        SAMRAI::hier::GlobalId patch;
        int dataId;
        std::vector<std::size_t> contributions; // indices in contributions_, in sum order
        std::size_t own;                        // index of the patch's own contribution
    };

public:
    // called for each transaction delivering data to a local patch, while the schedule is built;
    // returns the handle the transaction stores its contribution with
    std::size_t add(SAMRAI::hier::Box const& dst_node, SAMRAI::hier::Box const& src_node,
                    std::shared_ptr<SAMRAI::hier::BoxOverlap> const& overlap, int const dataId)
    {
        if (finalized_)
            throw std::runtime_error("BorderSums: contribution added after finalize");
        if (dst_node.isPeriodicImage())
            throw std::runtime_error("BorderSums: periodic image as destination");

        auto const typedOverlap = std::dynamic_pointer_cast<Overlap_t const>(overlap);
        if (!typedOverlap)
            throw std::runtime_error("BorderSums: unexpected overlap type");

        auto const size = detail::streamSize(*typedOverlap);
        // src_node is the periodic image of the contributing patch if the overlap is shifted
        contributions_.push_back(
            {phare_box_from<dimension>(src_node), typedOverlap, storeSize_, size});
        storeSize_ += size;

        destinationOf_(dst_node, dataId).contributions.push_back(contributions_.size() - 1);
        return contributions_.size() - 1;
    }


    // called once the schedule is built: adds each patch's own contribution, over the union of
    // the overlaps of the others, and sorts the contributions of each patch.
    // `boxOf(GlobalId)`: the box of a local patch
    void finalize(auto&& boxOf)
    {
        for (auto& destination : destinations_)
        {
            SAMRAI::hier::Box const& patchBox = boxOf(destination.patch);
            auto const own                    = ownOverlap_(destination, patchBox.getDim());
            auto const size                   = detail::streamSize(*own);
            contributions_.push_back({phare_box_from<dimension>(patchBox), own, storeSize_, size});
            storeSize_ += size;
            destination.own = contributions_.size() - 1;
            destination.contributions.push_back(destination.own);

            std::ranges::sort(destination.contributions, [&](auto const a, auto const b) {
                return core::lexicographicLess(contributions_[a].source, contributions_[b].source);
            });

            // one contribution per contributing patch or periodic image: an equal key would add
            // the same deposit twice
            auto const sameSource = [&](auto const a, auto const b) {
                return contributions_[a].source == contributions_[b].source;
            };
            auto const& sorted = destination.contributions;
            if (std::ranges::adjacent_find(sorted, sameSource) != sorted.end())
                throw std::runtime_error("BorderSums: two contributions from one patch box");
        }
        store_.resize(storeSize_);
        finalized_ = true;
    }

    // where the transaction of `handle` stores its contribution, in stream layout
    std::span<value_type> contribution(std::size_t const handle)
    {
        auto const& c = contributions_[handle];
        return {store_.data() + c.offset, c.size};
    }


    // after the schedule has run: the sum on each local patch, in place.
    // `dataOf(GlobalId, dataId)`: the FieldData_t of a local patch
    void fold(auto&& dataOf)
    {
        PHARE_LOG_SCOPE(3, "BorderSums::fold");

        for (auto const& destination : destinations_)
        {
            FieldData_t& data = dataOf(destination.patch, destination.dataId);
            auto const& own   = contributions_[destination.own];

            data.packTo(contribution(destination.own), *own.overlap);
            data.template unpackFrom<detail::SetZero<value_type>>(contribution(destination.own),
                                                                  *own.overlap);
            for (auto const i : destination.contributions)
                data.template unpackFrom<core::PlusEquals<value_type>>(contribution(i),
                                                                       *contributions_[i].overlap);
        }
    }

private:
    Destination& destinationOf_(SAMRAI::hier::Box const& dst_node, int const dataId)
    {
        auto const key = std::make_pair(dst_node.getGlobalId(), dataId);
        if (auto const it = destinationIndex_.find(key); it != destinationIndex_.end())
            return destinations_[it->second];
        destinationIndex_.emplace(key, destinations_.size());
        return destinations_.emplace_back(Destination{dst_node.getGlobalId(), dataId, {}, 0});
    }

    // the union of the overlaps of the other contributions, as disjoint boxes, per component
    std::shared_ptr<Overlap_t const> ownOverlap_(Destination const& destination,
                                                 SAMRAI::tbox::Dimension const& dim) const
    {
        using Components = decltype(detail::overlapComponents(std::declval<Overlap_t const&>()));
        auto constexpr N = std::tuple_size_v<Components>;

        std::array<SAMRAI::hier::BoxContainer, N> unions;
        for (auto const i : destination.contributions)
        {
            auto const components = detail::overlapComponents(*contributions_[i].overlap);
            for (std::size_t c = 0; c < N; ++c)
                for (auto const& box : components[c]->getDestinationBoxContainer())
                {
                    SAMRAI::hier::BoxContainer newPart{box};
                    newPart.removeIntersections(unions[c]);
                    unions[c].spliceBack(newPart);
                }
        }
        return detail::makeOverlap<Overlap_t>(
            unions, SAMRAI::hier::Transformation{SAMRAI::hier::IntVector::getZero(dim)});
    }

    std::vector<Contribution> contributions_;
    std::vector<Destination> destinations_;
    std::map<std::pair<SAMRAI::hier::GlobalId, int>, std::size_t> destinationIndex_;
    std::vector<value_type> store_;
    std::size_t storeSize_ = 0;
    bool finalized_        = false;
};



// BorderSums::finalize and fold on the patches of `level`
template<typename FieldData_t>
void finalizeBorderSums(BorderSums<FieldData_t>& sums, SAMRAI::hier::PatchLevel const& level)
{
    sums.finalize([&](auto const& id) -> auto const& { return level.getPatch(id)->getBox(); });
}

template<typename FieldData_t>
void foldBorderSums(BorderSums<FieldData_t>& sums, SAMRAI::hier::PatchLevel const& level)
{
    sums.fold([&](auto const& id, int const dataId) -> FieldData_t& {
        return dynamic_cast<FieldData_t&>(*level.getPatch(id)->getPatchData(dataId));
    });
}



/** A border-sum transaction that stores the contribution it carries in BorderSums instead of
 * adding it: on reception for a remote source, by copy from the source patch for a local one.
 * As a source, it packs like FieldBorderOpTransaction.
 */
template<typename FieldData_t>
class FieldBorderSumTransaction
    : public FieldBorderOpTransaction<FieldData_t,
                                      core::PlusEquals<typename FieldData_t::value_type>>
{
    using Super
        = FieldBorderOpTransaction<FieldData_t, core::PlusEquals<typename FieldData_t::value_type>>;

public:
    FieldBorderSumTransaction(std::shared_ptr<SAMRAI::hier::PatchLevel> const& dst_level,
                              std::shared_ptr<SAMRAI::hier::PatchLevel> const& src_level,
                              std::shared_ptr<SAMRAI::hier::BoxOverlap> const& overlap,
                              SAMRAI::hier::Box const& dst_node, SAMRAI::hier::Box const& src_node,
                              SAMRAI::xfer::RefineClasses::Data const** refine_data,
                              int const item_id, std::shared_ptr<BorderSums<FieldData_t>> sums,
                              std::size_t const handle)
        : Super{dst_level, src_level, overlap, dst_node, src_node, refine_data, item_id}
        , sums_{std::move(sums)}
        , handle_{handle}
    {
    }

    void unpackStream(SAMRAI::tbox::MessageStream& stream) override
    {
        PHARE_LOG_SCOPE(3, "FieldBorderSumTransaction::unpackStream");
        auto const contribution = sums_->contribution(handle_);
        stream.unpack(contribution.data(), contribution.size());
    }

    void copyLocalData() override
    {
        PHARE_LOG_SCOPE(3, "FieldBorderSumTransaction::copyLocalData");
        auto const& source = dynamic_cast<FieldData_t const&>(
            *this->d_src_level->getPatch(this->d_src_node.getGlobalId())
                 ->getPatchData(this->d_refine_data[this->d_item_id]->d_src));
        source.packTo(sums_->contribution(handle_), *this->d_overlap);
    }

private:
    std::shared_ptr<BorderSums<FieldData_t>> sums_;
    std::size_t handle_;
};



template<typename FieldData_t>
class FieldBorderSumTransactionFactory : public SAMRAI::xfer::RefineTransactionFactory
{
public:
    explicit FieldBorderSumTransactionFactory(std::shared_ptr<BorderSums<FieldData_t>> sums)
        : sums_{std::move(sums)}
    {
    }

    std::shared_ptr<SAMRAI::tbox::Transaction>
    allocate(std::shared_ptr<SAMRAI::hier::PatchLevel> const& dst_level,
             std::shared_ptr<SAMRAI::hier::PatchLevel> const& src_level,
             std::shared_ptr<SAMRAI::hier::BoxOverlap> const& overlap,
             SAMRAI::hier::Box const& dst_node, SAMRAI::hier::Box const& src_node,
             SAMRAI::xfer::RefineClasses::Data const** refine_data, int item_id,
             SAMRAI::hier::Box const& /*box*/, bool /*use_time_interpolation*/) const override
    {
        PHARE_LOG_SCOPE(3, "FieldBorderSumTransactionFactory::allocate");

        // only transactions delivering to a local patch carry a contribution to store
        bool const toLocalPatch
            = dst_node.getOwnerRank() == dst_level->getBoxLevel()->getMPI().getRank();
        auto const handle = toLocalPatch ? sums_->add(dst_node, src_node, overlap,
                                                      refine_data[item_id]->d_scratch)
                                         : std::size_t{0};

        return std::make_shared<FieldBorderSumTransaction<FieldData_t>>(
            dst_level, src_level, overlap, dst_node, src_node, refine_data, item_id, sums_, handle);
    }

    void
    preprocessScratchSpace(std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*level*/,
                           double /*fill_time*/,
                           SAMRAI::hier::ComponentSelector const& /*preprocess*/) const override
    {
    }

private:
    std::shared_ptr<BorderSums<FieldData_t>> sums_;
};


} // namespace PHARE::amr

#endif
