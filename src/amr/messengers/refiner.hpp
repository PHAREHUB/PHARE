#ifndef PHARE_REFINER_HPP
#define PHARE_REFINER_HPP

#include "communicator.hpp"
#include "core/data/vecfield/vecfield.hpp"

#include "amr/messengers/field_sum_transaction.hpp"

#include <tuple>
#include <cstdint>
#include <stdexcept>


namespace PHARE::amr
{

enum class RefinerType : std::uint16_t {
    GhostField = 0,
    PatchGhostField,
    InitField,
    InitInteriorPart,
    LevelBorderParticles,
    PatchFieldBorderSum,
    ExteriorGhostParticles,
    COUNT
};


auto constexpr static inline as_signed(RefinerType rt)
{
    return static_cast<int>(rt);
}
auto constexpr static inline refiner_type_enum_count()
{
    return static_cast<std::size_t>(RefinerType::COUNT);
}

template<typename ResourcesManager, RefinerType Type>
class Refiner : private Communicator<RefinerTypes, ResourcesManager::dimension>
{
    using FieldData_t = typename ResourcesManager::UserField_t::patch_data_type;


public:
    void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                       std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
    {
        auto levelNumber = level->getLevelNumber();

        for (auto& algo : this->algos)
        {
            // for GhostField we need schedules that take on the level where there is an
            // overlap (there is always for patches lying inside the level) and goes to
            // coarser level where there is not (patch lying on the level border) Create a
            // communication schedule that communicates data within a single level and
            // interpolates data from coarser hierarchy levels where needed.

            // Data will be communicated from the interiors of the source data on the given
            // level to the interiors and ghosts of destination data on the same level where
            // those sources and destinations overlap. Where they do not overlap,
            // data will be interpolated from source data on coarser levels in the patch
            // hierarchy. Data is time interpolated between old and new sources on coarser
            // levels when and where time interpolation is needed and copied from the source
            // components on the patch level into the destination components otherwise. Note
            // that the next coarser level number must correspond to a level in the
            // hierarchy that represents a region of coarser index space than the
            // destination level. Note that the schedule remains valid as long as the levels
            // involved in its creation do not change; thus, it can be used for multiple
            // data communication cycles.
            if constexpr (Type == RefinerType::GhostField)
            {
                this->add(algo,
                          algo->createSchedule(level, level->getNextCoarserHierarchyLevelNumber(),
                                               hierarchy),
                          levelNumber);
            }

            // the following schedule will only fill patch ghost nodes
            // not level border ghosts
            else if constexpr (Type == RefinerType::PatchGhostField)
            {
                this->add(algo, algo->createSchedule(level), levelNumber);
            }

            // schedule used to += density and flux for populations
            // on incomplete overlaped ghost box nodes
            else if constexpr (Type == RefinerType::PatchFieldBorderSum)
            {
                this->add(algo,
                          algo->createSchedule(
                              level, 0,
                              std::make_shared<FieldBorderSumTransactionFactory<FieldData_t>>()),
                          levelNumber);
            }

            // this createSchedule overload is used to initialize fields.
            // note that here we must take that createsSchedule() overload and put nullptr
            // as src since we want to take from coarser level everywhere. using the
            // createSchedule overload that takes level, next_coarser_level only would
            // result in interior ghost nodes to be filled with interior of neighbor patches
            // but there is nothing there.
            else if constexpr (Type == RefinerType::InitField)
            {
                this->add(algo, algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy),
                          levelNumber);
            }


            // here we create the schedule that will intialize the particles that lie within
            // the interior of the patches (no ghost, no coarse to fine). We take almost the
            // same overload as for fields above but the version that takes a
            // PatchLevelFillPattern. Here the PatchLevelInteriorFillPattern is used because
            // we want to fill particles only within the interior of the patches of the
            // level. The reason is that filling the their ghost regions with refined
            // particles would not ensure the ghosts to be clones of neighbor patches
            // particles if the splitting from coarser levels is not deterministic.
            else if constexpr (Type == RefinerType::InitInteriorPart)
            {
                this->add(algo,
                          algo->createSchedule(
                              std::make_shared<SAMRAI::xfer::PatchLevelInteriorFillPattern>(),
                              level, nullptr, levelNumber - 1, hierarchy),
                          levelNumber);
            }

            // here we create a schedule that will refine particles from coarser level and
            // put them into the level coarse to fine boundary. These are the
            // levelGhostParticlesOld particles. we thus take the same createSchedule
            // overload as above but pass it a PatchLevelBorderFillPattern.
            else if constexpr (Type == RefinerType::LevelBorderParticles)
            {
                this->add(algo,
                          algo->createSchedule(
                              std::make_shared<SAMRAI::xfer::PatchLevelBorderFillPattern>(), level,
                              nullptr, levelNumber - 1, hierarchy),
                          levelNumber);
            }


            else if constexpr (Type == RefinerType::ExteriorGhostParticles)
            {
                this->add(algo, algo->createSchedule(level), levelNumber);
            }
        }
    }



    void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                int const levelNumber, std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                double const initDataTime)
    {
        for (auto& algo : this->algos)
        {
            auto const& level = hierarchy->getPatchLevel(levelNumber);

            if constexpr (Type == RefinerType::InitInteriorPart)
            {
                auto schedule = algo->createSchedule(
                    std::make_shared<SAMRAI::xfer::PatchLevelInteriorFillPattern>(), level,
                    oldLevel, level->getNextCoarserHierarchyLevelNumber(), hierarchy);
                schedule->fillData(initDataTime);
            }
            else
            {
                auto schedule = algo->createSchedule(
                    level, oldLevel, level->getNextCoarserHierarchyLevelNumber(), hierarchy);
                schedule->fillData(initDataTime);
            }
        }
    }


    void fill(int const levelNumber, double const initDataTime) const
    {
        if (this->algos.size() == 0)
            throw std::runtime_error("Algorithms are not configured");

        for (auto const& algo : this->algos)
            this->findSchedule(algo, levelNumber)->fillData(initDataTime);
    }

    template<typename VecFieldT>
    void fill(VecFieldT& vec, int const levelNumber, double const fillTime)
    {
        for (auto const& algo : this->algos)
            this->findSchedule(algo, levelNumber)->fillData(fillTime);
    }


    /**
     * @Brief This overload creates a Refiner for communication with both spatial and
     * time interpolation. Data is communicated from the model vector field defined at
     * time t=n+1 and its version at time t=n (oldModel), onto the `ghost` vector field.
     *
     *
     * @param ghost represents the VecField that needs its ghost nodes filled
     * @param model represents the VecField from which data is taken (at
     * time t_coarse+dt_coarse)
     * @param oldModel represents the model VecField from which data is taken
     * at time t_coarse
     * @param rm is the ResourcesManager
     * @param refineOp is the spatial refinement operator
     * @param timeOp is the time interpolator
     *
     * @return the function returns a Refiner
     */
    Refiner(core::VecFieldNames const& ghost, core::VecFieldNames const& model,
            core::VecFieldNames const& oldModel, std::shared_ptr<ResourcesManager> const& rm,
            std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
            std ::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp,
            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> variableFillPattern = nullptr)
    {
        constexpr auto dimension = ResourcesManager::dimension;

        register_time_interpolated_vector_field( //
            rm, ghost, ghost, oldModel, model, refineOp, timeOp, variableFillPattern);
    }



    /**
     * @brief creates a Refiner for a scalar quantity with time refinement
     */
    Refiner(std::string const& ghost, std::string const& model, std::string const& oldModel,
            std::shared_ptr<ResourcesManager> const& rm,
            std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
            std ::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp,
            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> variableFillPattern = nullptr)
    {
        constexpr auto dimension = ResourcesManager::dimension;

        register_time_interpolated_resource( //
            rm, ghost, ghost, oldModel, model, refineOp, timeOp, variableFillPattern);
    }


    /**
     * @brief this overload creates a Refiner for communication without time interpolation
     * and from one quantity to the same quantity. It is typically used for initialization.
     */
    Refiner(core::VecFieldNames const& src_dest, std::shared_ptr<ResourcesManager> const& rm,
            std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
        : Refiner(src_dest, src_dest, rm, refineOp)
    {
    }


    /**
     * @brief this overload creates a Refiner for communication without time interpolation
     * and from one quantity to another quantity.
     */
    Refiner(core::VecFieldNames const& dst, core::VecFieldNames const& src,
            std::shared_ptr<ResourcesManager> const& rm,
            std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> variableFillPattern = nullptr)
    {
        register_vector_field(rm, dst, src, refineOp, variableFillPattern);
    }



    Refiner(std::string const& dst, std::string const& src,
            std::shared_ptr<ResourcesManager> const& rm,
            std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr)
    {
        auto&& [idDst, idSrc] = rm->getIDsList(dst, src);
        this->add_algorithm()->registerRefine(idDst, idSrc, idDst, refineOp, fillPattern);
    }



    /**
     * @brief This overload of makeRefiner creates a Refiner for communication from one
     * scalar quantity to itself without time interpolation.
     */
    Refiner(std::string const& name, std::shared_ptr<ResourcesManager> const& rm,
            std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
        : Refiner{name, name, rm, refineOp}
    {
    }



    auto& register_resource(auto& rm, auto& dst, auto& src, auto& scratch, auto&&... args)
    {
        auto&& [idDst, idSrc, idScrtch] = rm->getIDsList(dst, src, scratch);
        this->add_algorithm()->registerRefine(idDst, idSrc, idScrtch, args...);
        return *this;
    }


    auto& register_time_interpolated_resource(auto& rm, auto& dst, auto& src, auto& told,
                                              auto& tnew, auto&&... args)
    {
        auto&& [idDst, idSrc, idTold, idTnew] = rm->getIDsList(dst, src, told, tnew);
        this->add_algorithm()->registerRefine(idDst, idSrc, idTold, idTnew, idDst, args...);
        return *this;
    }


    auto& register_vector_field(auto& rm, auto& dst, auto& src, auto& refOp, auto& fillPat)
    {
        return (*this)
            .register_resource(rm, dst.xName, src.xName, dst.xName, refOp, fillPat)
            .register_resource(rm, dst.yName, src.yName, dst.yName, refOp, fillPat)
            .register_resource(rm, dst.zName, src.zName, dst.zName, refOp, fillPat);
    }


    auto& register_time_interpolated_vector_field(auto& rm, auto& dst, auto& src, auto& told,
                                                  auto& tnew, auto&&... args)
    {
        return (*this)
            .register_time_interpolated_resource(rm, dst.xName, src.xName, told.xName, tnew.xName,
                                                 args...)
            .register_time_interpolated_resource(rm, dst.yName, src.yName, told.yName, tnew.yName,
                                                 args...)
            .register_time_interpolated_resource(rm, dst.zName, src.zName, told.zName, tnew.zName,
                                                 args...);
    }
};



} // namespace PHARE::amr

#endif
