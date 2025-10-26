#ifndef PHARE_TEST_TAG_STRATEGY_HPP
#define PHARE_TEST_TAG_STRATEGY_HPP

#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include "amr/data/field/field_data.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"

#include <map>
#include <string>
#include <memory>


using namespace PHARE::core;
using namespace PHARE::amr;



/* When determining which DataFactory to use to create the correct geometry
 *  Samrai compares existing items within the RefinerPool
 *   See: SAMRAI::xfer::RefineClasses::getEquivalenceClassIndex"
 *   If we do not specify a separate VariableFillPattern for each item sent to
 *   Algo->registerRefine : we run the risk of all items only using the first registered
 *   DataFactory, thus all items will receive the Geometry of that item.
 *   We "hack" this as there is a typeid()== check on the variablefill pattern types
 */

class RhoVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class PVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};



class BXVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class BYVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class BZVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};


class EXVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class EYVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class EZVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};


class JXVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class JYVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class JZVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};




template<typename GridLayoutT, typename FieldT>
class TagStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy
{
public:
    TagStrategy(std::map<std::string, int> const& dataToAllocate,
                std::shared_ptr<SAMRAI::hier::RefineOperator>& refineOperator)
        : dataToAllocate_{dataToAllocate}
        , refineOp_{refineOperator}
    {
        auto registerRefine = [this](int id, auto& fillPattern) {
            algorithm_.registerRefine(id, id, id, refineOp_, fillPattern);
        };

        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> BxVariableFillPattern
            = std::make_shared<BXVariableFillPattern>();
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> ByVariableFillPattern
            = std::make_shared<BYVariableFillPattern>();
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> BzVariableFillPattern
            = std::make_shared<BZVariableFillPattern>();

        registerRefine(dataToAllocate_["Bx"], BxVariableFillPattern);
        registerRefine(dataToAllocate_["By"], ByVariableFillPattern);
        registerRefine(dataToAllocate_["Bz"], BzVariableFillPattern);

        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> ExVariableFillPattern
            = std::make_shared<EXVariableFillPattern>();
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> EyVariableFillPattern
            = std::make_shared<EYVariableFillPattern>();
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> EzVariableFillPattern
            = std::make_shared<EZVariableFillPattern>();

        registerRefine(dataToAllocate_["Ex"], ExVariableFillPattern);
        registerRefine(dataToAllocate_["Ey"], EyVariableFillPattern);
        registerRefine(dataToAllocate_["Ez"], EzVariableFillPattern);

        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> JxVariableFillPattern
            = std::make_shared<JXVariableFillPattern>();
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> JyVariableFillPattern
            = std::make_shared<JYVariableFillPattern>();
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> JzVariableFillPattern
            = std::make_shared<JZVariableFillPattern>();

        registerRefine(dataToAllocate_["Jx"], JxVariableFillPattern);
        registerRefine(dataToAllocate_["Jy"], JyVariableFillPattern);
        registerRefine(dataToAllocate_["Jz"], JzVariableFillPattern);

        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> rhoVariableFillPattern
            = std::make_shared<RhoVariableFillPattern>();

        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> pVariableFillPattern
            = std::make_shared<PVariableFillPattern>();

        registerRefine(dataToAllocate_["Rho"], rhoVariableFillPattern);
        registerRefine(dataToAllocate_["P"], pVariableFillPattern);
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

                    for (auto const bix : layout.AMRGhostBoxFor(field))
                    {
                        auto position  = layout.fieldNodeCoordinates(field, bix);
                        auto const lix = layout.AMRToLocal(bix);
                        field(lix)     = affineFill(position, dataId);
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

    static double affineFill(Point<double, GridLayoutT::dimension> position, int dataId)
    {
        static auto constexpr dim = GridLayoutT::dimension;

        // parameter for linear function ax + by + cz + d
        double a = dataId + 1.0;
        double b = dataId + 10.0;
        double c = dataId + 100.0;
        double d = dataId + 1000.0;

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
