#ifndef PHARE_AMR_MESSENGERS_SCHEDULER_HPP
#define PHARE_AMR_MESSENGERS_SCHEDULER_HPP

#include "refiner.hpp"


#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>

#include <SAMRAI/xfer/VariableFillPattern.h>


namespace PHARE::amr
{


struct Scheduler
{
};


struct RefinerScheduler : public Scheduler
{
    constexpr static std::size_t n_refiner_types = refiner_type_enum_count(); //
    using Schedule                               = SAMRAI::xfer::RefineSchedule;
    using Algorithm                              = SAMRAI::xfer::RefineAlgorithm;

    struct Duo
    {
        std::vector<std::unique_ptr<Algorithm>> a;
        std::map<int, std::vector<std::shared_ptr<Schedule>>> s;
        std::vector<std::function<void(SAMRAI::hier::PatchLevel&, double /*fillTime*/)>> funcs;
    };



    template<auto rtype, typename Extra = void>
    RefinerScheduler& add(std::string const& dst,
                          std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                          std::shared_ptr<SAMRAI::hier::PatchLevel> const& level);


    template<auto rtype>
    auto& add_algorithm(std::string const& dst,
                        std::function<void(SAMRAI::hier::PatchLevel&, double)> optional = {})
    {
        auto& duo = duos[dst][as_signed(rtype)];
        duo.a.emplace_back(std::make_unique<Algorithm>());

        if (bool{optional})
            duo.funcs.emplace_back(optional);

        else // else default to regular
            duo.funcs.emplace_back([=, this](auto& lvl, double time) {
                (*this).template call<rtype>(dst, lvl, time);
            });

        return *this;
    }



    template<auto rtype>
    auto& register_resource(auto& rm, auto& dst, auto& src, auto& scratch, auto&&... args)
    {
        auto&& [idDst, idSrc, idScrtch] = rm->getIDsList(dst, src, scratch);

        auto& duo = duos[dst][as_signed(rtype)];
        if (!duo.a.size())
            add_algorithm<rtype>(dst);
        duo.a.back()->registerRefine(idDst, idSrc, idScrtch, args...);
        return *this;
    }

    template<auto rtype>
    auto& register_time_interpolated_resource(auto& rm, auto& dst, auto& src, auto& told,
                                              auto& tnew, auto&&... args)
    {
        auto&& [idDst, idSrc, idTold, idTnew] = rm->getIDsList(dst, src, told, tnew);
        auto& duo                             = duos[dst][as_signed(rtype)];
        if (!duo.a.size())
            add_algorithm<rtype>(dst);
        duo.a.back()->registerRefine(idDst, idSrc, idTold, idTnew, idDst, args...);
        return *this;
    }

    template<auto rtype>
    auto& register_vector_field(auto& rm, auto& dst, auto& src, auto&& refOp, auto&& fillPat)
    {
        register_resource<rtype>(rm, dst.xName, src.xName, dst.xName, refOp, fillPat);
        register_resource<rtype>(rm, dst.yName, src.yName, dst.yName, refOp, fillPat);
        register_resource<rtype>(rm, dst.zName, src.zName, dst.zName, refOp, fillPat);
        return *this;
    }

    template<auto rtype>
    auto& register_tensor_field(auto& rm, auto& dst, auto& src, auto&& fillPat)
    {
        for (std::size_t i = 0; i < dst.componentNames().size(); ++i)
            register_resource<rtype>(rm, dst.componentNames()[i], src.componentNames()[i],
                                     dst.componentNames()[i], /*refOp*/ nullptr, fillPat);
        return *this;
    }

    template<auto rtype>
    auto& register_time_interpolated_vector_field(auto& rm, auto& dst, auto& src, auto& told,
                                                  auto& tnew, auto&&... args)
    {
        register_time_interpolated_resource<rtype>(rm, dst.xName, src.xName, told.xName, tnew.xName,
                                                   args...);
        register_time_interpolated_resource<rtype>(rm, dst.yName, src.yName, told.yName, tnew.yName,
                                                   args...);
        register_time_interpolated_resource<rtype>(rm, dst.zName, src.zName, told.zName, tnew.zName,
                                                   args...);
        return *this;
    }



    // fill functions are general purpose
    template<auto rtype>
    void fill(std::string const& dst, SAMRAI::hier::PatchLevel& level, double const initDataTime)
    {
        auto& duo = duos[dst][as_signed(rtype)];
        for (auto& func : duo.funcs)
            func(level, initDataTime);
    }

    template<auto rtype> // override for split schedules
    void fill(std::string const& dst, SAMRAI::hier::PatchLevel& level, double const initDataTime,
              int const idx)
    {
        auto& duo = duos[dst][as_signed(rtype)];
        duo.funcs[idx](level, initDataTime);
    }


    // call functions deal with schedules directly - for advanced usage!
    template<auto rtype>
    void call(std::string const& dst, int const levelNumber, double const initDataTime)
    {
        auto& duo = duos[dst][as_signed(rtype)];
        assert(duo.a.size() == duo.s[levelNumber].size());
        for (std::size_t i = 0; i < duo.a.size(); ++i)
            duo.s[levelNumber][i]->fillData(initDataTime);
    }

    template<auto rtype>
    void call(std::string const& dst, SAMRAI::hier::PatchLevel const& level,
              double const initDataTime)
    {
        call<rtype>(dst, level.getLevelNumber(), initDataTime);
    }

    template<auto rtype> // override for split schedules
    void call(std::string const& dst, SAMRAI::hier::PatchLevel const& level,
              double const initDataTime, int const idx)
    {
        duos[dst][as_signed(rtype)].s[level.getLevelNumber()][idx]->fillData(initDataTime);
    }


    std::unordered_map<std::string, std::array<Duo, n_refiner_types>> duos;
};


template<auto Type, typename Extra>
RefinerScheduler&
RefinerScheduler::add(std::string const& dst,
                      std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                      std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
{
    auto const levelNumber = level->getLevelNumber();

    auto& duo = duos[dst][as_signed(Type)];
    assert(duo.a.size());

    for (auto& algo : duo.a)
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
            duo.s[levelNumber].emplace_back(algo->createSchedule(
                level, level->getNextCoarserHierarchyLevelNumber(), hierarchy));
        }

        // the following schedule will only fill patch ghost nodes
        // not level border ghosts
        else if constexpr (Type == RefinerType::PatchGhostField)
        {
            duo.s[levelNumber].emplace_back(algo->createSchedule(level));
        }

        else if constexpr (Type == RefinerType::PatchFieldBorderSum)
        {
            static_assert(not std::is_same_v<Extra, void> && "NO!");

            duo.s[levelNumber].emplace_back(algo->createSchedule(
                level, 0, std::make_shared<FieldBorderSumTransactionFactory<Extra>>()));
        }

        // this createSchedule overload is used to initialize fields.
        // note that here we must take that createsSchedule() overload and put nullptr
        // as src since we want to take from coarser level everywhere. using the
        // createSchedule overload that takes level, next_coarser_level only would
        // result in interior ghost nodes to be filled with interior of neighbor patches
        // but there is nothing there.
        else if constexpr (Type == RefinerType::InitField)
        {
            duo.s[levelNumber].emplace_back(
                algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy));
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
            duo.s[levelNumber].emplace_back(algo->createSchedule(
                std::make_shared<SAMRAI::xfer::PatchLevelInteriorFillPattern>(), level, nullptr,
                levelNumber - 1, hierarchy));
        }

        // here we create a schedule that will refine particles from coarser level and
        // put them into the level coarse to fine boundary. These are the
        // levelGhostParticlesOld particles. we thus take the same createSchedule
        // overload as above but pass it a PatchLevelBorderFillPattern.
        else if constexpr (Type == RefinerType::LevelBorderParticles)
        {
            duo.s[levelNumber].emplace_back(
                algo->createSchedule(std::make_shared<SAMRAI::xfer::PatchLevelBorderFillPattern>(),
                                     level, nullptr, levelNumber - 1, hierarchy));
        }


        else if constexpr (Type == RefinerType::ExteriorGhostParticles)
        {
            duo.s[levelNumber].emplace_back(algo->createSchedule(level));
        }
    }

    return *this;
}



} // namespace PHARE::amr

#endif
