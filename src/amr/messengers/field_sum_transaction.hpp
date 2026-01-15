#ifndef PHARE_AMR_MESSENGERS_FIELD_SUM_TRANSACTION_HPP
#define PHARE_AMR_MESSENGERS_FIELD_SUM_TRANSACTION_HPP

#include "core/logger.hpp"

#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include <stdexcept>

namespace PHARE::amr
{


/** * @brief FieldBorderSumTransaction is used to += pop density and flux on ghost box overlaps
 *
 * A FieldBorderSumTransaction is a SAMRAI Transaction created by the
 * FieldBorderSumTransactionFactory provided (via createShedule) to schedules that accumulate
 * incomplete density and flux on ghost box overlaps.
 *
 * Due to the lack of neighbor particle contributions, some domain nodes and ghost nodes
 * have incomplete moments after deposition. The complement of these nodes is what has
 * been deposited on (also incomplete) neighbor nodes.
 *
 * Default SAMRAI transaction calls PatchData::copy and PatchData::packStream
 * This transaction defines these override to these methods to call specific methods
 * of FieldData to perform the += instead of =.
 * These methods are copyAndSum and unpackStreamAndSum.
 *
 */
template<typename FieldData_t>
class FieldBorderSumTransaction : public SAMRAI::tbox::Transaction
{
public:
    FieldBorderSumTransaction(std::shared_ptr<SAMRAI::hier::PatchLevel> const& dst_level,
                              std::shared_ptr<SAMRAI::hier::PatchLevel> const& src_level,
                              std::shared_ptr<SAMRAI::hier::BoxOverlap> const& overlap,
                              SAMRAI::hier::Box const& dst_node, SAMRAI::hier::Box const& src_node,
                              SAMRAI::xfer::RefineClasses::Data const** refine_data, int item_id)
        : d_dst_level(dst_level)
        , d_src_level(src_level)
        , d_overlap(overlap)
        , d_dst_node(dst_node)
        , d_src_node(src_node)
        , d_refine_data(refine_data)
        , d_item_id(item_id)
        , d_incoming_bytes(0)
        , d_outgoing_bytes(0)
    {
        TBOX_ASSERT(dst_level);
        TBOX_ASSERT(src_level);
        TBOX_ASSERT(overlap);
        TBOX_ASSERT(dst_node.getLocalId() >= 0);
        TBOX_ASSERT(src_node.getLocalId() >= 0);
        TBOX_ASSERT(item_id >= 0);
        TBOX_ASSERT(refine_data[item_id] != 0);

        TBOX_ASSERT_OBJDIM_EQUALITY4(*dst_level, *src_level, dst_node, src_node);
    }

    virtual ~FieldBorderSumTransaction() {}


    virtual bool canEstimateIncomingMessageSize();

    virtual size_t computeIncomingMessageSize();

    virtual size_t computeOutgoingMessageSize();

    virtual int getSourceProcessor();

    virtual int getDestinationProcessor();

    virtual void packStream(SAMRAI::tbox::MessageStream& stream);

    virtual void unpackStream(SAMRAI::tbox::MessageStream& stream);

    virtual void copyLocalData();

    virtual void printClassData(std::ostream& stream) const;

private:
    std::shared_ptr<SAMRAI::hier::PatchLevel> d_dst_level;
    std::shared_ptr<SAMRAI::hier::PatchLevel> d_src_level;
    std::shared_ptr<SAMRAI::hier::BoxOverlap> d_overlap;
    SAMRAI::hier::Box d_dst_node;
    SAMRAI::hier::Box d_src_node;
    SAMRAI::xfer::RefineClasses::Data const** d_refine_data;
    int d_item_id;
    size_t d_incoming_bytes;
    size_t d_outgoing_bytes;
};


template<typename FieldData_t>
bool FieldBorderSumTransaction<FieldData_t>::canEstimateIncomingMessageSize()
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::canEstimateIncomingMessageSize");
    bool can_estimate = false;
    if (getSourceProcessor() == d_src_level->getBoxLevel()->getMPI().getRank())
    {
        can_estimate = d_src_level->getPatch(d_src_node.getGlobalId())
                           ->getPatchData(d_refine_data[d_item_id]->d_src)
                           ->canEstimateStreamSizeFromBox();
    }
    else
    {
        can_estimate = d_dst_level->getPatch(d_dst_node.getGlobalId())
                           ->getPatchData(d_refine_data[d_item_id]->d_scratch)
                           ->canEstimateStreamSizeFromBox();
    }
    return can_estimate;
}


template<typename FieldData_t>
size_t FieldBorderSumTransaction<FieldData_t>::computeIncomingMessageSize()
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::computeIncomingMessageSize");
    d_incoming_bytes = d_dst_level->getPatch(d_dst_node.getGlobalId())
                           ->getPatchData(d_refine_data[d_item_id]->d_scratch)
                           ->getDataStreamSize(*d_overlap);
    return d_incoming_bytes;
}

template<typename FieldData_t>
size_t FieldBorderSumTransaction<FieldData_t>::computeOutgoingMessageSize()
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::computeOutgoingMessageSize");
    d_outgoing_bytes = d_src_level->getPatch(d_src_node.getGlobalId())
                           ->getPatchData(d_refine_data[d_item_id]->d_src)
                           ->getDataStreamSize(*d_overlap);
    return d_outgoing_bytes;
}

template<typename FieldData_t>
int FieldBorderSumTransaction<FieldData_t>::getSourceProcessor()
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::getSourceProcessor");
    return d_src_node.getOwnerRank();
}

template<typename FieldData_t>
int FieldBorderSumTransaction<FieldData_t>::getDestinationProcessor()
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::getDestinationProcessor");
    return d_dst_node.getOwnerRank();
}

template<typename FieldData_t>
void FieldBorderSumTransaction<FieldData_t>::packStream(SAMRAI::tbox::MessageStream& stream)
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::packStream");
    d_src_level->getPatch(d_src_node.getGlobalId())
        ->getPatchData(d_refine_data[d_item_id]->d_src)
        ->packStream(stream, *d_overlap);
}

template<typename FieldData_t>
void FieldBorderSumTransaction<FieldData_t>::unpackStream(SAMRAI::tbox::MessageStream& stream)
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::unpackStream");
    std::shared_ptr<FieldData_t> onode_dst_data(
        SAMRAI_SHARED_PTR_CAST<FieldData_t, SAMRAI::hier::PatchData>(
            d_dst_level->getPatch(d_dst_node.getGlobalId())
                ->getPatchData(d_refine_data[d_item_id]->d_scratch)));
    TBOX_ASSERT(onode_dst_data);

    onode_dst_data->unpackStreamAndSum(stream, *d_overlap);
}


template<typename FieldData_t>
void FieldBorderSumTransaction<FieldData_t>::printClassData(std::ostream& stream) const
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::printClassData");
    throw std::runtime_error("FieldBorderSumTransaction::printClassData!");
}

template<typename FieldData_t>
void FieldBorderSumTransaction<FieldData_t>::copyLocalData()
{
    PHARE_LOG_SCOPE(2, "FieldBorderSumTransaction::copyLocalData");
    std::shared_ptr<FieldData_t> onode_dst_data(
        SAMRAI_SHARED_PTR_CAST<FieldData_t, SAMRAI::hier::PatchData>(
            d_dst_level->getPatch(d_dst_node.getGlobalId())
                ->getPatchData(d_refine_data[d_item_id]->d_scratch)));
    TBOX_ASSERT(onode_dst_data);

    std::shared_ptr<FieldData_t> onode_src_data(
        SAMRAI_SHARED_PTR_CAST<FieldData_t, SAMRAI::hier::PatchData>(
            d_src_level->getPatch(d_src_node.getGlobalId())
                ->getPatchData(d_refine_data[d_item_id]->d_src)));
    TBOX_ASSERT(onode_src_data);

    onode_dst_data->sum(*onode_src_data, *d_overlap);
}


template<typename FieldData_t>
class FieldBorderSumTransactionFactory : public SAMRAI::xfer::RefineTransactionFactory
{
public:
    std::shared_ptr<SAMRAI::tbox::Transaction>
    allocate(std::shared_ptr<SAMRAI::hier::PatchLevel> const& dst_level,
             std::shared_ptr<SAMRAI::hier::PatchLevel> const& src_level,
             std::shared_ptr<SAMRAI::hier::BoxOverlap> const& overlap,
             SAMRAI::hier::Box const& dst_node, SAMRAI::hier::Box const& src_node,
             SAMRAI::xfer::RefineClasses::Data const** refine_data, int item_id,
             SAMRAI::hier::Box const& box, bool use_time_interpolation) const override
    {
        NULL_USE(box);
        NULL_USE(use_time_interpolation);

        TBOX_ASSERT(dst_level);
        TBOX_ASSERT(src_level);
        TBOX_ASSERT(overlap);
        TBOX_ASSERT(dst_node.getLocalId() >= 0);
        TBOX_ASSERT(src_node.getLocalId() >= 0);
        TBOX_ASSERT(refine_data != 0);
        TBOX_ASSERT_OBJDIM_EQUALITY4(*dst_level, *src_level, dst_node, src_node);

        PHARE_LOG_SCOPE(2, "FieldBorderSumTransactionFactory::allocate");
        return std::make_shared<FieldBorderSumTransaction<FieldData_t>>(
            dst_level, src_level, overlap, dst_node, src_node, refine_data, item_id);
    }

    void
    preprocessScratchSpace(std::shared_ptr<SAMRAI::hier::PatchLevel> const& level, double fill_time,
                           SAMRAI::hier::ComponentSelector const& preprocess_vector) const override
    {
        PHARE_LOG_SCOPE(2, "FieldBorderSumTransactionFactory::preprocessScratchSpace");

        // noop
    }
};

} // namespace PHARE::amr

#endif // PHARE_AMR_MESSENGERS_FIELD_SUM_TRANSACTION_HPP
