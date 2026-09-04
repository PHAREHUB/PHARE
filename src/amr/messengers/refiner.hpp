#ifndef PHARE_REFINER_HPP
#define PHARE_REFINER_HPP

#include "communicator.hpp"

#include "amr/messengers/field_operate_transaction.hpp"
#include "core/utilities/types.hpp"

#include <stdexcept>


namespace PHARE::amr
{

enum class RefinerType {
    GhostField,
    InitField,
    InitInteriorPart,
    LevelBorderField,
    LevelBorderParticles,
    PatchGhostField,
    PatchIonPopBorderSum,
    PatchIonPopBorderMax,
    PatchTensorFieldBorderSum,
    PatchVecFieldBorderMax,
    ExteriorGhostParticles
};



template<typename ResourcesManager, typename Model, RefinerType Type>
class Refiner : private Communicator<RefinerTypes, ResourcesManager::dimension>
{
    using UserTypes_t = ResourcesManager::template ResourceUserTypesFor<Model>;

    using FieldData_t = UserTypes_t::UserField_t::patch_data_type;

    using SetMaxOp     = core::SetMax<typename FieldData_t::value_type>;
    using PlusEqualsOp = core::PlusEquals<typename FieldData_t::value_type>;

    // hard coded rank cause there's no real tensorfields that use this code yet
    using TensorFieldData_t = UserTypes_t::template UserTensorField_t<2>::patch_data_type;
    using VecFieldData_t    = UserTypes_t::template UserTensorField_t<1>::patch_data_type;

    std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat_ = nullptr;

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
                                               hierarchy, patchStrat_.get()),
                          levelNumber);
            }

            else if constexpr (Type == RefinerType::PatchGhostField)
            {
                this->add(algo, algo->createSchedule(level, patchStrat_.get()), levelNumber);
            }


            // schedule used to += flux and density for populations on incomplete overlapped
            // ghost box nodes. Flux (VecFieldData_t) and both density fields (FieldData_t) can
            // be registered together onto this single algorithm (see add_resource()); the
            // transaction factory resolves each item's concrete PatchData type at runtime so
            // they are all communicated together in a single set of MPI messages.
            else if constexpr (Type == RefinerType::PatchIonPopBorderSum)
            {
                this->add(
                    algo,
                    algo->createSchedule(level, patchStrat_.get(),
                                         std::make_shared<MultiFieldBorderOpTransactionFactory<
                                             PlusEqualsOp, FieldData_t, VecFieldData_t>>()),
                    levelNumber);
            }


            else if constexpr (Type == RefinerType::PatchTensorFieldBorderSum)
            {
                this->add(
                    algo,
                    algo->createSchedule(
                        level, patchStrat_.get(),
                        std::make_shared<
                            FieldBorderOpTransactionFactory<TensorFieldData_t, PlusEqualsOp>>()),
                    levelNumber);
            }


            // schedule used to == max of density and bulk velocity on complete overlapped
            // ghost box nodes. Density (FieldData_t) and bulk velocity (VecFieldData_t) can be
            // registered together onto this single algorithm (see add_resource()); the
            // transaction factory resolves each item's concrete PatchData type at runtime so
            // they are all communicated together in a single set of MPI messages.
            else if constexpr (Type == RefinerType::PatchIonPopBorderMax)
            {
                this->add(
                    algo,
                    algo->createSchedule(
                        level, patchStrat_.get(),
                        std::make_shared<MultiFieldBorderOpTransactionFactory<SetMaxOp, FieldData_t,
                                                                              VecFieldData_t>>()),
                    levelNumber);
            }


            // schedule used to == max of density and flux for populations
            // on complete overlaped ghost box nodes
            else if constexpr (Type == RefinerType::PatchVecFieldBorderMax)
            {
                this->add(algo,
                          algo->createSchedule(
                              level, patchStrat_.get(),
                              std::make_shared<
                                  FieldBorderOpTransactionFactory<VecFieldData_t, SetMaxOp>>()),
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
                this->add(algo,
                          algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy,
                                               patchStrat_.get()),
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
                              level, nullptr, levelNumber - 1, hierarchy, patchStrat_.get()),
                          levelNumber);
            }


            // this will fill ghosts on level border only (i.e. not patch ghost nodes)
            // it is not used at the moment but would be useful if the moment ghosts
            // are filled on level border (see Note in hybhyb)
            else if constexpr (Type == RefinerType::LevelBorderField)
            {
                this->add(algo,
                          algo->createSchedule(
                              std::make_shared<SAMRAI::xfer::PatchLevelBorderFillPattern>(), level,
                              level->getNextCoarserHierarchyLevelNumber(), hierarchy,
                              patchStrat_.get()),
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
                              nullptr, levelNumber - 1, hierarchy, patchStrat_.get()),
                          levelNumber);
            }


            else if constexpr (Type == RefinerType::ExteriorGhostParticles)
            {
                this->add(algo, algo->createSchedule(level, patchStrat_.get()), levelNumber);
            }

            else
            {
                throw std::runtime_error("No Schedule for RefinerType");
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
                    oldLevel, level->getNextCoarserHierarchyLevelNumber(), hierarchy,
                    patchStrat_.get());
                schedule->fillData(initDataTime);
            }
            else
            {
                auto schedule = algo->createSchedule(level, oldLevel,
                                                     level->getNextCoarserHierarchyLevelNumber(),
                                                     hierarchy, patchStrat_.get());
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
     * @brief creates a Refiner for a scalar quantity with time refinement
     */
    Refiner(std::string const& ghost, std::string const& model, std::string const& oldModel,
            std::shared_ptr<ResourcesManager> const& rm,
            std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
            std ::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp,
            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> variableFillPattern = nullptr,
            std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat          = nullptr)
    {
        constexpr auto dimension = ResourcesManager::dimension;

        patchStrat_ = patchStrat;

        register_time_interpolated_resource( //
            rm, ghost, ghost, oldModel, model, refineOp, timeOp, variableFillPattern);
    }




    Refiner(std::string const& dst, std::string const& src,
            std::shared_ptr<ResourcesManager> const& rm,
            std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr,
            std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat  = nullptr)
    {
        patchStrat_ = patchStrat;

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


    /**
     * @brief creates an empty Refiner with no resource registered yet. Resources are added one
     * at a time via add_resource(), all sharing the same underlying algorithm/schedule so that
     * they end up being communicated together in a single set of MPI messages.
     */
    Refiner() = default;


    /**
     * @brief adds a resource to the single shared algorithm of this Refiner, alongside any
     * resource previously added via add_resource(). Unlike register_resource(), which creates a
     * new algorithm (and thus a new schedule) on every call, this reuses the same algorithm so
     * that heterogeneous resources (e.g. a scalar Field and a VecField) fuse into one schedule.
     */
    auto& add_resource(auto& rm, auto& dst, auto& src, auto& scratch, auto&&... args)
    {
        if (this->algos.empty())
            this->add_algorithm();
        auto&& [idDst, idSrc, idScrtch] = rm->getIDsList(dst, src, scratch);
        this->algos[0]->registerRefine(idDst, idSrc, idScrtch, args...);
        return *this;
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
