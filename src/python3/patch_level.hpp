#ifndef PHARE_PYTHON_PATCH_LEVEL_HPP
#define PHARE_PYTHON_PATCH_LEVEL_HPP



#include "amr/wrappers/hierarchy.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"

#include "python3/patch_data.hpp"



#include <string>
#include <cstring>
#include <cstddef>



namespace PHARE::pydata
{


template<typename Model, typename Enable = void>
class PatchLevel;

template<typename Model_t>
class AnyPatchLevel
{
    using GridLayout = Model_t::gridlayout_type;

    static constexpr std::size_t dimension = GridLayout::dimension;

public:
    AnyPatchLevel(amr::Hierarchy& hierarchy, Model_t& model, std::size_t lvl)
        : lvl(lvl)
        , hierarchy{hierarchy}
        , model{model}
    {
    }


protected:
    auto getField(auto& field)
    {
        std::vector<PatchData<py_array_t<double>, dimension>> patchDatas;
        auto visit = [&](auto& grid, auto patchID, auto /*ilvl*/) {
            setPyPatchDataFromField(patchDatas.emplace_back(), field, grid, patchID);
        };
        amr::visitLevel<GridLayout>(*hierarchy.getPatchLevel(lvl), rm, visit, field);
        return patchDatas;
    }
    auto getVecFieldComponent(auto& vecfield, auto const component)
    {
        std::vector<PatchData<py_array_t<double>, dimension>> patchDatas;
        auto visit = [&](auto& grid, auto patchID, auto /*ilvl*/) {
            setPyPatchDataFromField(patchDatas.emplace_back(), vecfield(component), grid, patchID);
        };
        amr::visitLevel<GridLayout>(*hierarchy.getPatchLevel(lvl), rm, visit, vecfield);
        return patchDatas;
    }

    auto static vecfield_component(std::string const& component)
    {
        return core::Components::componentMap().at(component);
    }

protected:
    std::size_t lvl;
    amr::Hierarchy& hierarchy;
    Model_t& model;
    Model_t::resources_manager_type& rm = *model.resourcesManager;
};


template<typename Model>
class __attribute__((visibility("hidden")))
PatchLevel<Model, std::enable_if_t<solver::is_hybrid_model_v<Model>>> : AnyPatchLevel<Model>
{
    using Super = AnyPatchLevel<Model>;

public:
    using Model_t                          = Model;
    using GridLayout                       = Model_t::gridlayout_type;
    static constexpr std::size_t dimension = GridLayout::dimension;

    PatchLevel(amr::Hierarchy& hierarchy, Model& model, std::size_t lvl)
        : Super(hierarchy, model, lvl)
    {
    }

    auto getB(std::string const& component)
    {
        return this->getVecFieldComponent(this->model.state.electromag.B,
                                          this->vecfield_component(component));
    }

    auto getE(std::string const& component)
    {
        return this->getVecFieldComponent(this->model.state.electromag.E,
                                          this->vecfield_component(component));
    }

    auto& getIonPop(auto const& popName)
    {
        for (auto& pop : this->model.state.ions)
            if (popName == pop.name())
                return pop;
        throw std::runtime_error("no population found named: " + popName);
    }

    auto getNi() { return this->getField(this->model.state.ions.massDensity()); }

    auto getN(std::string const& popName)
    {
        return this->getField(getIonPop(popName).chargeDensity());
    }

    auto getVi(std::string const& component)
    {
        return this->getVecFieldComponent(this->model.state.ions.velocity(),
                                          this->vecfield_component(component));
    }

    auto getFlux(std::string const& component, std::string const& popName)
    {
        return this->getVecFieldComponent(getIonPop(popName).flux(),
                                          this->vecfield_component(component));
    }

    auto getParticles(std::string const& popName)
    {
        std::vector<PatchData<core::ParticleArray<dimension>*, dimension>> patchDatas;

        auto& pop = getIonPop(popName);

        auto visit = [&](auto& grid, auto patchID, auto /*ilvl*/) {
            setPatchDataFromGrid(patchDatas.emplace_back(&pop.domainParticles()), grid, patchID);
        };

        amr::visitLevel<GridLayout>( //
            *this->hierarchy.getPatchLevel(this->lvl), this->rm, visit, pop);

        return patchDatas;
    }
};


template<typename Model>
class __attribute__((visibility("hidden")))
PatchLevel<Model, std::enable_if_t<solver::is_mhd_model_v<Model>>> : AnyPatchLevel<Model>
{
    using Super = AnyPatchLevel<Model>;

public:
    using Model_t    = Model;
    using GridLayout = Model_t::gridlayout_type;

    PatchLevel(amr::Hierarchy& hierarchy, Model& model, std::size_t const lvl)
        : Super(hierarchy, model, lvl)
    {
    }
};



} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PATCH_LEVEL_H*/
