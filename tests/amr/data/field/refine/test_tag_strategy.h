#ifndef PHARE_TEST_TAG_STRATEGY_H
#define PHARE_TEST_TAG_STRATEGY_H


#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include "amr/data/field/field_data.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/utilities/constants.h"
#include "core/utilities/point/point.h"


#include <map>
#include <string>

using namespace PHARE::core;
using namespace PHARE::amr;


template<typename GridLayoutT, typename FieldT>
class TagStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy
{
public:
    TagStrategy(std::map<std::string, int> const& dataToAllocate,
                std::shared_ptr<SAMRAI::hier::RefineOperator>& refineOperator)
        : dataToAllocate_{dataToAllocate}
        , refineOp_{refineOperator}
    {
        for (auto const& nameToIds : dataToAllocate_)
        {
            algorithm_.registerRefine(nameToIds.second, nameToIds.second, nameToIds.second,
                                      refineOp_);
        }
    }

    virtual ~TagStrategy() = default;

    void initializeLevelData(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             int const levelNumber, double const initDataTime, bool const,
                             bool const,
                             std::shared_ptr<SAMRAI::hier::PatchLevel> const& = std::shared_ptr<
                                 SAMRAI::hier::PatchLevel>(),
                             bool const allocateData = true) override
    {
        if (allocateData)
        {
            auto level = hierarchy->getPatchLevel(levelNumber);
            for (auto& patch : *level)
            {
                for (auto const& dataPair : dataToAllocate_)
                {
                    auto const& dataId = dataPair.second;
                    patch->allocatePatchData(dataId, initDataTime);
                }
            }
        }

        if (levelNumber == 0)
        {
            auto level = hierarchy->getPatchLevel(levelNumber);
            for (auto& patch : *level)
            {
                for (auto const& variablesId : dataToAllocate_)
                {
                    auto const& dataId = variablesId.second;
                    auto fieldData     = std::dynamic_pointer_cast<FieldData<GridLayoutT, FieldT>>(
                        patch->getPatchData(dataId));

                    auto& layout = fieldData->gridLayout;
                    auto& field  = fieldData->field;

                    if constexpr (dim == 1)
                    {
                        std::uint32_t gsi_X = layout.ghostStartIndex(field, Direction::X);
                        std::uint32_t gei_X = layout.ghostEndIndex(field, Direction::X);

                        for (std::uint32_t ix = gsi_X; ix <= gei_X; ++ix)
                        {
                            auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix);
                            field(ix)     = affineFill(position);
                        }
                    }
                    if constexpr (dim == 2)
                    {
                        std::uint32_t gsi_X = layout.ghostStartIndex(field, Direction::X);
                        std::uint32_t gei_X = layout.ghostEndIndex(field, Direction::X);
                        std::uint32_t gsi_Y = layout.ghostStartIndex(field, Direction::Y);
                        std::uint32_t gei_Y = layout.ghostEndIndex(field, Direction::Y);

                        for (std::uint32_t ix = gsi_X; ix <= gei_X; ++ix)
                        {
                            for (std::uint32_t iy = gsi_Y; iy <= gei_Y; ++iy)
                            {
                                auto position
                                    = layout.fieldNodeCoordinates(field, layout.origin(), ix, iy);
                                auto affineVal = affineFill(position);
                                field(ix, iy)  = affineFill(position);
                            }
                        }
                    }
                    if constexpr (dim == 3)
                    {
                        std::uint32_t gsi_X = layout.ghostStartIndex(field, Direction::X);
                        std::uint32_t gei_X = layout.ghostEndIndex(field, Direction::X);
                        std::uint32_t gsi_Y = layout.ghostStartIndex(field, Direction::Y);
                        std::uint32_t gei_Y = layout.ghostEndIndex(field, Direction::Y);
                        std::uint32_t gsi_Z = layout.ghostStartIndex(field, Direction::Z);
                        std::uint32_t gei_Z = layout.ghostEndIndex(field, Direction::Z);

                        for (std::uint32_t ix = gsi_X; ix <= gei_X; ++ix)
                        {
                            for (std::uint32_t iy = gsi_Y; iy <= gei_Y; ++iy)
                            {
                                for (std::uint32_t iz = gsi_Z; iz <= gei_Z; ++iz)
                                {
                                    auto position = layout.fieldNodeCoordinates(
                                        field, layout.origin(), ix, iy, iz);
                                    field(ix, iy, iz) = affineFill(position);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            /* // create schedule */
            auto refineSchedule = algorithm_.createSchedule(hierarchy->getPatchLevel(levelNumber),
                                                            nullptr, levelNumber - 1, hierarchy);

            refineSchedule->fillData(0.);
        }
    }

    void resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const&,
                                     int const, int const) override
    {
        // do nothing
    }

    static double affineFill(Point<double, GridLayoutT::dimension> position)
    {
        static auto constexpr dim = GridLayoutT::dimension;

        // parameter for linear function ax + by + cz + d
        double a = 1.0;
        double b = 10.0;
        double c = 100.0;
        double d = 1000.0;

        if constexpr (dim == 1)
        {
            return a * position[dirX] + d;
        }
        if constexpr (dim == 2)
        {
            return a * position[dirX] + b * position[dirY] + d;
        }
        if constexpr (dim == 3)
        {
            return a * position[dirX] + b * position[dirY] + c * position[dirZ] + d;
        }
    }



private:
    std::map<std::string, int> dataToAllocate_;
    std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp_;
    SAMRAI::xfer::RefineAlgorithm algorithm_;
    static auto constexpr dim = GridLayoutT::dimension;
};

#endif
