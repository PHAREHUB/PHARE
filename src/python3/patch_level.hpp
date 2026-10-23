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


template<typename Model>
class PatchLevel;

template<typename Model_t>
class AnyPatchLevel
{
    static constexpr std::size_t dimension = Model_t::dimension;

public:
    AnyPatchLevel(amr::Hierarchy& hierarchy, Model_t& model, std::size_t lvl)
        : lvl(lvl)
        , hierarchy{hierarchy}
        , model{model}
    {
    }


protected:
    template<typename GridLayout>
    auto getField(auto& field)
    {
        return getPatchDatasForLevel<GridLayout>(hierarchy, model, lvl, field);
    }

    template<typename GridLayout>
    auto getVecFieldComponent(auto& vecfield, auto const component)
    {
        return getPatchDatasForLevel<GridLayout>(hierarchy, model, lvl, vecfield,
                                                 [=](auto& vf) -> auto& { return vf(component); });
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
    requires(solver::is_hybrid_model_v<Model>)
class __attribute__((visibility("hidden"))) PatchLevel<Model> : AnyPatchLevel<Model>
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
        return this->template getVecFieldComponent<GridLayout>(this->model.state.electromag.B,
                                                               this->vecfield_component(component));
    }

    auto getE(std::string const& component)
    {
        return this->template getVecFieldComponent<GridLayout>(this->model.state.electromag.E,
                                                               this->vecfield_component(component));
    }

    auto& getIonPop(auto const& popName)
    {
        for (auto& pop : this->model.state.ions)
            if (popName == pop.name())
                return pop;
        throw std::runtime_error("no population found named: " + popName);
    }

    auto getNi()
    {
        return this->template getField<GridLayout>(this->model.state.ions.massDensity());
    }

    auto getN(std::string const& popName)
    {
        return this->template getField<GridLayout>(getIonPop(popName).chargeDensity());
    }

    auto getVi(std::string const& component)
    {
        return this->template getVecFieldComponent<GridLayout>(this->model.state.ions.velocity(),
                                                               this->vecfield_component(component));
    }

    auto getFlux(std::string const& component, std::string const& popName)
    {
        return this->template getVecFieldComponent<GridLayout>(getIonPop(popName).flux(),
                                                               this->vecfield_component(component));
    }

    auto getParticles(std::string const& popName)
    {
        using ParticleArray_t = Model_t::particle_array_type;
        std::vector<PatchData<ParticleArray_t*, dimension>> patchDatas;

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
    requires(solver::is_mhd_model_v<Model>)
class __attribute__((visibility("hidden"))) PatchLevel<Model> : AnyPatchLevel<Model>
{
    using Super = AnyPatchLevel<Model>;

public:
    using Model_t    = Model;
    using GridLayout = Model_t::gridlayout_type;

    PatchLevel(amr::Hierarchy& hierarchy, Model& model, std::size_t const lvl)
        : Super(hierarchy, model, lvl)
    {
    }

    auto getRho() { return this->template getField<GridLayout>(this->model.state.rho); }
    auto getP() { return this->template getField<GridLayout>(this->model.state.P); }
    auto getEtot() { return this->template getField<GridLayout>(this->model.state.Etot); }

    auto getV(std::string const& component)
    {
        return this->template getVecFieldComponent<GridLayout>(this->model.state.V,
                                                               this->vecfield_component(component));
    }
    auto getB(std::string const& component)
    {
        return this->template getVecFieldComponent<GridLayout>(this->model.state.B,
                                                               this->vecfield_component(component));
    }
    auto getRhoV(std::string const& component)
    {
        return this->template getVecFieldComponent<GridLayout>(this->model.state.rhoV,
                                                               this->vecfield_component(component));
    }
    auto getE(std::string const& component)
    {
        return this->template getVecFieldComponent<GridLayout>(this->model.state.E,
                                                               this->vecfield_component(component));
    }
    auto getJ(std::string const& component)
    {
        return this->template getVecFieldComponent<GridLayout>(this->model.state.J,
                                                               this->vecfield_component(component));
    }
};



} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PATCH_LEVEL_H*/
