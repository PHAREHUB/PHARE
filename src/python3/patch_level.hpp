#ifndef PHARE_PYTHON_PATCH_LEVEL_HPP
#define PHARE_PYTHON_PATCH_LEVEL_HPP

#include "phare_solver.hpp"
#include "python3/patch_data.hpp"

#include <string>
#include <cstring>
#include <cstddef>

#include "python3/patch_data.hpp"


namespace PHARE::pydata
{
template<auto opts>
class __attribute__((visibility("hidden"))) PatchLevel
{
public:
    static constexpr std::size_t dimension = opts.dimension;

    using PHARESolverTypes = solver::PHARE_Types<opts>;
    using HybridModel      = PHARESolverTypes::HybridModel_t;
    using GridLayout       = HybridModel::gridlayout_type;

    PatchLevel(amr::Hierarchy& hierarchy, HybridModel& model, std::size_t lvl)
        : lvl_(lvl)
        , hierarchy_{hierarchy}
        , model_{model}
    {
    }

    auto getDensity()
    {
        std::vector<PatchData<std::vector<double>, dimension>> patchDatas;
        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            setPatchDataFromField(patchDatas.emplace_back(), ions.chargeDensity(), grid, patchID);
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, ions);

        return patchDatas;
    }

    auto getPopDensities()
    {
        using Inner = decltype(getDensity());

        std::unordered_map<std::string, Inner> pop_data;
        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            for (auto const& pop : ions)
            {
                if (!pop_data.count(pop.name()))
                    pop_data.emplace(pop.name(), Inner());

                setPatchDataFromField(pop_data.at(pop.name()).emplace_back(), pop.chargeDensity(),
                                      grid, patchID);
            }
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, ions);

        return pop_data;
    }

    template<typename VecField, typename Map>
    auto getVecFields(VecField& vecField, Map& container, GridLayout& grid, std::string patchID,
                      std::string outer_key)
    {
        if (!container.count(outer_key))
            container.emplace(outer_key, typename Map::mapped_type());

        auto& inner = container.at(outer_key);

        for (auto& [id, type] : core::Components::componentMap())
        {
            auto& field = vecField.getComponent(type);

            if (!inner.count(field.name()))
                inner.emplace(field.name(), decltype(getDensity())());

            setPatchDataFromField(inner.at(field.name()).emplace_back(), field, grid, patchID);
        }
    }

    auto getBulkVelocity()
    {
        decltype(getPopDensities()) bulkV;

        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            for (auto& [id, type] : core::Components::componentMap())
            {
                auto& field = ions.velocity().getComponent(type);

                if (!bulkV.count(field.name()))
                    bulkV.emplace(field.name(), decltype(getDensity())());

                setPatchDataFromField(bulkV.at(field.name()).emplace_back(), field, grid, patchID);
            }
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, ions);

        return bulkV;
    }

    auto getPopFlux()
    {
        using Inner = decltype(getPopDensities());

        std::unordered_map<std::string, Inner> pop_data;

        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            for (auto const& pop : ions)
                getVecFields(pop.flux(), pop_data, grid, patchID, pop.name());
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, ions);

        return pop_data;
    }

    auto getEM()
    {
        using Inner = decltype(getPopDensities());

        std::unordered_map<std::string, Inner> em_data;

        auto& em = model_.state.electromag;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            for (auto& vecFieldPtr : {&em.B, &em.E})
            {
                getVecFields(*vecFieldPtr, em_data, grid, patchID, vecFieldPtr->name());
            }
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, em);

        return em_data;
    }


    auto getB(std::string componentName)
    {
        std::vector<PatchData<std::vector<double>, dimension>> patchDatas;

        auto& B = model_.state.electromag.B;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            auto compo = PHARE::core::Components::componentMap().at(componentName);
            setPatchDataFromField(patchDatas.emplace_back(), B.getComponent(compo), grid, patchID);
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, B);

        return patchDatas;
    }


    auto getE(std::string componentName)
    {
        std::vector<PatchData<std::vector<double>, dimension>> patchDatas;

        auto& E = model_.state.electromag.E;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            auto compo = PHARE::core::Components::componentMap().at(componentName);
            setPatchDataFromField(patchDatas.emplace_back(), E.getComponent(compo), grid, patchID);
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, E);

        return patchDatas;
    }



    auto getVi(std::string componentName)
    {
        std::vector<PatchData<std::vector<double>, dimension>> patchDatas;

        auto& V = model_.state.ions.velocity();

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            auto compo = PHARE::core::Components::componentMap().at(componentName);
            setPatchDataFromField(patchDatas.emplace_back(), V.getComponent(compo), grid, patchID);
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, V);

        return patchDatas;
    }



    auto getPopFluxCompo(std::string component, std::string popName)
    {
        std::vector<PatchData<std::vector<double>, dimension>> patchDatas;

        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            auto compo = PHARE::core::Components::componentMap().at(component);
            for (auto const& pop : ions)
                if (pop.name() == popName)
                    setPatchDataFromField(patchDatas.emplace_back(), pop.flux().getComponent(compo),
                                          grid, patchID);
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, ions);

        return patchDatas;
    }




    auto getEx() { return getE("x"); }
    auto getEy() { return getE("y"); }
    auto getEz() { return getE("z"); }

    auto getBx() { return getB("x"); }
    auto getBy() { return getB("y"); }
    auto getBz() { return getB("z"); }

    auto getVix() { return getVi("x"); }
    auto getViy() { return getVi("y"); }
    auto getViz() { return getVi("z"); }

    auto getFx(std::string pop) { return getPopFluxCompo("x", pop); }
    auto getFy(std::string pop) { return getPopFluxCompo("y", pop); }
    auto getFz(std::string pop) { return getPopFluxCompo("z", pop); }




    auto getParticles(std::string userPopName)
    {
        using Nested = std::vector<PatchData<core::ContiguousParticles<dimension>, dimension>>;
        using Inner  = std::unordered_map<std::string, Nested>;

        std::unordered_map<std::string, Inner> pop_particles;

        auto getParticleData = [&](Inner& inner, GridLayout& grid, std::string patchID,
                                   std::string key, auto& particles) {
            if (particles.size() == 0)
                return;

            if (!inner.count(key))
                inner.emplace(key, Nested());

            auto& patch_data = inner[key].emplace_back(particles.size());
            setPatchDataFromGrid(patch_data, grid, patchID);
            core::ParticlePacker<dimension>{particles}.pack(patch_data.data);
        };

        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            for (auto& pop : ions)
            {
                if ((userPopName != "" and userPopName == pop.name()) or userPopName == "all")
                {
                    if (!pop_particles.count(pop.name()))
                        pop_particles.emplace(pop.name(), Inner());

                    auto& inner = pop_particles.at(pop.name());

                    getParticleData(inner, grid, patchID, "domain", pop.domainParticles());
                    getParticleData(inner, grid, patchID, "levelGhost", pop.levelGhostParticles());
                }
            }
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, ions);

        return pop_particles;
    }

private:
    std::size_t lvl_;
    amr::Hierarchy& hierarchy_;
    HybridModel& model_;
};


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PATCH_LEVEL_H*/
