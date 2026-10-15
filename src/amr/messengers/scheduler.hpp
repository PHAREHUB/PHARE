#ifndef PHARE_AMR_MESSENGERS_SCHEDULER_HPP
#define PHARE_AMR_MESSENGERS_SCHEDULER_HPP

#include "amr/messengers/refiner.hpp"


#include <map>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <functional>
#include <unordered_map>

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

    using PatchStrategy = SAMRAI::xfer::RefinePatchStrategy;

    struct Duo
    {
        std::vector<std::unique_ptr<Algorithm>> a;
        std::vector<PatchStrategy*> ps; // per algorithm, non owning, may be null
        std::map<int, std::vector<std::shared_ptr<Schedule>>> s;
        std::vector<std::function<void(SAMRAI::hier::PatchLevel const&, double /*fillTime*/)>>
            funcs;
    };



    // (re)creates the schedules of all algorithms of `dst` for `level`
    template<auto rtype, typename ResourcesManager>
    RefinerScheduler& add(std::string const& dst,
                          std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                          std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
    {
        auto& duo = duos[dst][as_signed(rtype)];
        assert(duo.a.size());

        auto& schedules = duo.s[level->getLevelNumber()];
        schedules.clear(); // level (re)registration replaces previous schedules
        for (std::size_t i = 0; i < duo.a.size(); ++i)
            schedules.emplace_back(make_refine_schedule<rtype, ResourcesManager>(
                *duo.a[i], hierarchy, level, duo.ps[i]));

        return *this;
    }


    // fills the new level `levelNumber` from `oldLevel` and coarser levels, with temporary
    // schedules
    template<auto rtype>
    void regrid(std::string const& dst,
                std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                int const levelNumber, std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                double const initDataTime)
    {
        auto& duo = duos[dst][as_signed(rtype)];
        for (std::size_t i = 0; i < duo.a.size(); ++i)
            make_regrid_schedule<rtype>(*duo.a[i], hierarchy, levelNumber, oldLevel, duo.ps[i])
                ->fillData(initDataTime);
    }


    template<auto rtype>
    auto& add_algorithm(std::string const& dst,
                        std::function<void(SAMRAI::hier::PatchLevel const&, double)> optional = {},
                        PatchStrategy* patchStrat = nullptr)
    {
        auto& duo      = duos[dst][as_signed(rtype)];
        auto const idx = duo.a.size();
        duo.a.emplace_back(std::make_unique<Algorithm>());
        duo.ps.emplace_back(patchStrat);

        if (bool{optional})
            duo.funcs.emplace_back(optional);

        else // else default to filling only the schedule of this algorithm
            duo.funcs.emplace_back([=, this](auto const& lvl, double time) {
                (*this).template call<rtype>(dst, lvl, time, idx);
            });

        return *this;
    }



    template<auto rtype>
    std::size_t n_algorithms(std::string const& dst)
    {
        return duos[dst][as_signed(rtype)].a.size();
    }



    // registers onto the last algorithm added for `key`, resources registered together
    // are communicated together by the same schedule
    template<auto rtype>
    auto& register_resource(std::string const& key, auto& rm, auto const& dst, auto const& src,
                            auto const& scratch, auto&&... args)
    {
        auto&& [idDst, idSrc, idScrtch] = rm->getIDsList(dst, src, scratch);

        auto& duo = duos[key][as_signed(rtype)];
        if (!duo.a.size())
            add_algorithm<rtype>(key);
        duo.a.back()->registerRefine(idDst, idSrc, idScrtch, args...);
        return *this;
    }

    template<auto rtype>
    auto& register_time_interpolated_resource(std::string const& key, auto& rm, auto const& dst,
                                              auto const& src, auto const& told,
                                              auto const& tnew, auto&&... args)
    {
        auto&& [idDst, idSrc, idTold, idTnew] = rm->getIDsList(dst, src, told, tnew);

        auto& duo = duos[key][as_signed(rtype)];
        if (!duo.a.size())
            add_algorithm<rtype>(key);
        duo.a.back()->registerRefine(idDst, idSrc, idTold, idTnew, idDst, args...);
        return *this;
    }




    // fill functions are general purpose
    template<auto rtype>
    void fill(std::string const& dst, SAMRAI::hier::PatchLevel const& level,
              double const initDataTime)
    {
        auto& duo = duos[dst][as_signed(rtype)];
        for (auto& func : duo.funcs)
            func(level, initDataTime);
    }

    template<auto rtype> // override for split schedules
    void fill(std::string const& dst, SAMRAI::hier::PatchLevel const& level,
              double const initDataTime, int const idx)
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


} // namespace PHARE::amr

#endif
