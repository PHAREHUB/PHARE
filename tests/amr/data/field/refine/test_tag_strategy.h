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
                        std::uint32_t iStartX = layout.ghostStartIndex(field, Direction::X);
                        std::uint32_t iEndX   = layout.ghostEndIndex(field, Direction::X);

                        for (std::uint32_t ix = iStartX; ix <= iEndX; ++ix)
                        {
                            auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix);
                            field(ix)     = affineFill(position);
                        }
                    }
                    if constexpr (dim == 2) {}
                    if constexpr (dim == 3) {}
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
        double a = 0.2;
        double b = 0.4;
        double c = 0.6;
        double d = 0.8;

        if constexpr (dim == 1)
        {
            return a * position[dirX] + b;
        }
        if constexpr (dim == 2) {}
        if constexpr (dim == 3) {}
    }



private:
    std::map<std::string, int> dataToAllocate_;
    std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp_;
    SAMRAI::xfer::RefineAlgorithm algorithm_;
    static auto constexpr dim = GridLayoutT::dimension;
};

#endif
