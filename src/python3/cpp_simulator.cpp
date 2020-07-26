#include "phare/include.h"

#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/chrono.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"

#include "phare/include.h"
#include "core/utilities/types.h"
#include "core/utilities/mpi_utils.h"
#include "core/data/particles/particle_packer.h"


namespace py = pybind11;

template<typename T>
using py_array_t = py::array_t<T, py::array::c_style | py::array::forcecast>;


namespace PHARE::pydata
{
template<typename Data, std::size_t dim>
struct PatchData
{
    static auto constexpr dimension = dim;
    std::string patchID;
    std::string origin;
    py_array_t<std::size_t> lower{dim};
    py_array_t<std::size_t> upper{dim};
    std::size_t nGhosts;
    Data data;

    PatchData() {}

    template<typename... Args>
    PatchData(Args&&... args)
        : data{std::forward<Args...>(args...)}
    {
    }
};


template<typename PatchData>
void setPatchData(PatchData& data, std::string patchID, std::string origin,
                  std::array<std::size_t, PatchData::dimension> lower,
                  std::array<std::size_t, PatchData::dimension> upper)
{
    constexpr std::size_t bytes = PatchData::dimension * sizeof(size_t);
    std::memcpy(data.lower.request().ptr, lower.data(), bytes);
    std::memcpy(data.upper.request().ptr, upper.data(), bytes);
    data.patchID = patchID;
    data.origin  = origin;
}

template<typename PatchData, typename GridLayout>
void setPatchDataFromGrid(PatchData& pdata, GridLayout& grid, std::string patchID)
{
    setPatchData(pdata, patchID, grid.origin().str(),
                 grid.AMRBox().lower.template toArray<std::size_t>(),
                 grid.AMRBox().upper.template toArray<std::size_t>());
}

template<typename PatchData, typename Field, typename GridLayout>
void setPatchDataFromField(PatchData& pdata, Field const& field, GridLayout& grid,
                           std::string patchID)
{
    setPatchDataFromGrid(pdata, grid, patchID);
    pdata.nGhosts = static_cast<std::size_t>(
        GridLayout::nbrGhosts(GridLayout::centering(field.physicalQuantity())[0]));
    pdata.data.assign(field.data(), field.data() + field.size());
}



template<std::size_t dim, std::size_t interpOrder, std::size_t nbrRefPart>
class PatchLevel
{
public:
    static constexpr std::size_t dimension     = dim;
    static constexpr std::size_t interp_order  = interpOrder;
    static constexpr std::size_t nbRefinedPart = nbrRefPart;

    using PHARETypes  = PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using HybridModel = typename PHARETypes::HybridModel_t;
    using GridLayout  = typename HybridModel::gridlayout_type;

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
            setPatchDataFromField(patchDatas.emplace_back(), ions.density(), grid, patchID);
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

                setPatchDataFromField(pop_data.at(pop.name()).emplace_back(), pop.density(), grid,
                                      patchID);
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

        for (auto& [id, type] : core::Components::componentMap)
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
            for (auto& [id, type] : core::Components::componentMap)
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
            auto compo = PHARE::core::Components::componentMap.at(componentName);
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
            auto compo = PHARE::core::Components::componentMap.at(componentName);
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
            auto compo = PHARE::core::Components::componentMap.at(componentName);
            setPatchDataFromField(patchDatas.emplace_back(), V.getComponent(compo), grid, patchID);
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, V);

        return patchDatas;
    }



    auto getPopFluxCompo(std::string component, std::string popName)
    {
        using Inner = decltype(getPopDensities());

        std::vector<PatchData<std::vector<double>, dimension>> patchDatas;

        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t /*iLevel*/) {
            auto compo = PHARE::core::Components::componentMap.at(component);
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
                    getParticleData(inner, grid, patchID, "patchGhost", pop.patchGhostParticles());
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

template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
class DataWrangler
{
public:
    using This                                 = DataWrangler;
    static constexpr std::size_t dimension     = _dimension;
    static constexpr std::size_t interp_order  = _interp_order;
    static constexpr std::size_t nbRefinedPart = _nbRefinedPart;

    using PHARETypes  = PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using HybridModel = typename PHARETypes::HybridModel_t;

    DataWrangler(std::shared_ptr<ISimulator> const& simulator,
                 std::shared_ptr<amr::Hierarchy> const& hierarchy)
        : simulator_{simulator}
        , hierarchy_{hierarchy}
    {
        auto simDict = initializer::PHAREDictHandler::INSTANCE().dict()["simulation"];
        if (!core::makeAtRuntime<Maker>(
                simDict["dimension"].template to<int>(), simDict["interp_order"].template to<int>(),
                simDict["refined_particle_nbr"].template to<int>(), Maker{*this}))
            throw std::runtime_error("Runtime diagnostic deduction failed");
    }


    auto getNumberOfLevels() const { return hierarchy_->getNumberOfLevels(); }

    auto getPatchLevel(size_t lvl)
    {
        return PatchLevel<_dimension, _interp_order, _nbRefinedPart>{
            *hierarchy_, *simulator_ptr_->getHybridModel(), lvl};
    }

    auto sort_merge_1d(std::vector<PatchData<std::vector<double>, dimension>> const&& input,
                       bool shared_patch_border = false)
    {
        std::vector<std::pair<double, const PatchData<std::vector<double>, dimension>*>> sorted;
        for (auto const& data : input)
            sorted.emplace_back(core::Point<double, 1>::fromString(data.origin)[0], &data);
        std::sort(sorted.begin(), sorted.end(), [](auto& a, auto& b) { return a.first < b.first; });
        std::vector<double> ret;
        for (size_t i = 0; i < sorted.size(); i++)
        { // skip empty patches in case of unequal patches across MPI domains
            if (!sorted[i].second->data.size())
                continue;
            auto& data   = sorted[i].second->data;
            auto& ghosts = sorted[i].second->nGhosts;
            auto end     = ghosts;
            // primal nodes share a cell wall when patches touch so drop duplicate value if so
            if (shared_patch_border)
                end = i == sorted.size() - 1 ? end : end + 1;
            ret.insert(std::end(ret), std::begin(data) + ghosts, std::end(data) - end);
        }
        return ret;
    }

    auto sync(std::vector<PatchData<std::vector<double>, dimension>> const& input)
    {
        int mpi_size = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        std::vector<PatchData<std::vector<double>, dimension>> collected;

        auto reinterpret_array = [](py_array_t<std::size_t> const& py_array) {
            return *reinterpret_cast<std::array<std::size_t, dimension>*>(py_array.request().ptr);
        };

        auto collect = [&](PatchData<std::vector<double>, dimension> const& patch_data) {
            auto patchIDs = core::mpi::collect(patch_data.patchID, mpi_size);
            auto origins  = core::mpi::collect(patch_data.origin, mpi_size);
            auto lower    = core::mpi::collect(reinterpret_array(patch_data.lower), mpi_size);
            auto upper    = core::mpi::collect(reinterpret_array(patch_data.upper), mpi_size);
            auto ghosts   = core::mpi::collect(patch_data.nGhosts, mpi_size);
            auto datas    = core::mpi::collect(patch_data.data, mpi_size);

            for (int i = 0; i < mpi_size; i++)
            {
                PatchData<std::vector<double>, dimension>& data = collected.emplace_back();
                setPatchData(data, patchIDs[i], origins[i], lower[i], upper[i]);
                data.nGhosts = ghosts[i];
                data.data    = std::move(datas[i]);
            }
        };

        std::size_t max = core::mpi::max(input.size(), mpi_size);

        PatchData<std::vector<double>, dimension> empty;

        for (size_t i = 0; i < max; i++)
        {
            if (i < input.size())
                collect(input[i]);
            else
                collect(empty);
        }
        return collected;
    }

    auto sync_merge(std::vector<PatchData<std::vector<double>, dimension>> const& input,
                    bool primal)
    {
        if constexpr (dimension == 1)
            return sort_merge_1d(sync(input), primal);

        static_assert("Unhandled dimension: sort_merge_fluid");
    }

private:
    std::shared_ptr<ISimulator> simulator_;
    std::shared_ptr<amr::Hierarchy> hierarchy_;
    Simulator<dimension, interp_order, nbRefinedPart>* simulator_ptr_ = nullptr;

    struct Maker
    {
        Maker(DataWrangler& _dw)
            : dw{_dw}
        {
        }

        template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
        bool operator()(std::size_t userDim, std::size_t userInterpOrder,
                        std::size_t userNbRefinedPart, Dimension dimension_fn,
                        InterpOrder interp_order_fn, NbRefinedPart nbRefinedPart_fn)
        {
            if (userDim == dimension_fn() and userInterpOrder == interp_order_fn()
                and userNbRefinedPart == nbRefinedPart_fn())
            {
                std::size_t constexpr d  = dimension_fn();
                std::size_t constexpr io = interp_order_fn();
                std::size_t constexpr nb = nbRefinedPart_fn();

                // extra if constexpr as cast is templated and not generic interface
                if constexpr (d == dimension and io == interp_order and nb == nbRefinedPart)
                    return (dw.simulator_ptr_
                            = dynamic_cast<PHARE::Simulator<d, io, nb>*>(dw.simulator_.get()));
            }
            return 0;
        }

        DataWrangler& dw;
    };
};



template<typename Type, std::size_t dimension>
void declarePatchData(py::module& m, std::string key)
{
    using PatchDataType = PatchData<Type, dimension>;
    py::class_<PatchDataType>(m, key.c_str())
        .def_readonly("patchID", &PatchDataType::patchID)
        .def_readonly("origin", &PatchDataType::origin)
        .def_readonly("lower", &PatchDataType::lower)
        .def_readonly("upper", &PatchDataType::upper)
        .def_readonly("nGhosts", &PatchDataType::nGhosts)
        .def_readonly("data", &PatchDataType::data);
}

template<std::size_t dim>
void declareDim(py::module& m)
{
    using CP         = core::ContiguousParticles<dim>;
    std::string name = "ContiguousParticles_" + std::to_string(dim);
    py::class_<CP, std::shared_ptr<CP>>(m, name.c_str())
        .def(py::init<std::size_t>())
        .def_readwrite("iCell", &CP::iCell)
        .def_readwrite("delta", &CP::delta)
        .def_readwrite("weight", &CP::weight)
        .def_readwrite("charge", &CP::charge)
        .def_readwrite("v", &CP::v)
        .def("size", &CP::size);

    name = "PatchData" + name;
    declarePatchData<CP, dim>(m, name.c_str());
}

template<typename Simulator, typename PyClass>
void declareSimulator(PyClass&& sim)
{
    sim.def("initialize", &Simulator::initialize)
        .def("advance", &Simulator::advance)
        .def("startTime", &Simulator::startTime)
        .def("currentTime", &Simulator::currentTime)
        .def("endTime", &Simulator::endTime)
        .def("timeStep", &Simulator::timeStep)
        .def("to_str", &Simulator::to_str)
        .def("domain_box", &Simulator::domainBox)
        .def("cell_width", &Simulator::cellWidth);
}


template<typename _dim, typename _interp, typename _nbRefinedPart>
void declare(py::module& m)
{
    constexpr auto dim           = _dim{}();
    constexpr auto interp        = _interp{}();
    constexpr auto nbRefinedPart = _nbRefinedPart{}();


    std::string type_string = "_" + std::to_string(dim) + "_" + std::to_string(interp) + "_"
                              + std::to_string(nbRefinedPart);

    using Sim        = Simulator<dim, interp, nbRefinedPart>;
    std::string name = "Simulator" + type_string;
    declareSimulator<Sim>(
        py::class_<Sim, std::shared_ptr<Sim>>(m, name.c_str())
            .def_property_readonly_static("dims", [](py::object) { return Sim::dimension; })
            .def_property_readonly_static("interp_order",
                                          [](py::object) { return Sim::interp_order; })
            .def_property_readonly_static("refined_particle_nbr",
                                          [](py::object) { return Sim::nbRefinedPart; }));


    name = "make_simulator" + type_string;
    m.def(name.c_str(), [](std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
        return std::shared_ptr<Sim>{std::move(makeSimulator<dim, interp, nbRefinedPart>(hier))};
    });
    m.def("make_diagnostic_manager",
          [](std::shared_ptr<Sim> const& sim, std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
              return std::make_shared<RuntimeDiagnosticInterface>(*sim, *hier);
          });

    using DW = DataWrangler<dim, interp, nbRefinedPart>;
    name     = "DataWrangler" + type_string;
    py::class_<DW, std::shared_ptr<DW>>(m, name.c_str())
        .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def(py::init<std::shared_ptr<ISimulator> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def("sync_merge", &DW::sync_merge)
        .def("getPatchLevel", &DW::getPatchLevel)
        .def("getNumberOfLevels", &DW::getNumberOfLevels);

    using PL = PatchLevel<dim, interp, nbRefinedPart>;
    name     = "PatchLevel_" + type_string;

    py::class_<PL, std::shared_ptr<PL>>(m, name.c_str())
        .def("getEM", &PL::getEM)
        .def("getE", &PL::getE)
        .def("getB", &PL::getB)
        .def("getBx", &PL::getBx)
        .def("getBy", &PL::getBy)
        .def("getBz", &PL::getBz)
        .def("getEx", &PL::getEx)
        .def("getEy", &PL::getEy)
        .def("getEz", &PL::getEz)
        .def("getVix", &PL::getVix)
        .def("getViy", &PL::getViy)
        .def("getViz", &PL::getViz)
        .def("getDensity", &PL::getDensity)
        .def("getBulkVelocity", &PL::getBulkVelocity)
        .def("getPopDensities", &PL::getPopDensities)
        .def("getPopFluxes", &PL::getPopFlux)
        .def("getFx", &PL::getFx)
        .def("getFy", &PL::getFy)
        .def("getFz", &PL::getFz)
        .def("getParticles", &PL::getParticles, py::arg("userPopName") = "all");

    using _Splitter
        = PHARE::amr::Splitter<_dim, _interp, PHARE::core::RefinedParticlesConst<nbRefinedPart>>;
    name = "Splitter" + type_string;
    py::class_<_Splitter, std::shared_ptr<_Splitter>>(m, name.c_str())
        .def(py::init<>())
        .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
        .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });
}


template<typename Dimension, typename InterpOrder, typename... NbRefinedParts>
void declare(py::module& m, std::tuple<Dimension, InterpOrder, NbRefinedParts...> const&)
{
    core::apply(std::tuple<NbRefinedParts...>{}, [&](auto& nbRefinedPart) {
        declare<Dimension, InterpOrder, std::decay_t<decltype(nbRefinedPart)>>(m);
    });
}


PYBIND11_MODULE(cpp, m)
{
    py::class_<SamraiLifeCycle, std::shared_ptr<SamraiLifeCycle>>(m, "SamraiLifeCycle")
        .def(py::init<>())
        .def("reset", &SamraiLifeCycle::reset);

    py::class_<PHARE::amr::Hierarchy, std::shared_ptr<PHARE::amr::Hierarchy>>(m, "AMRHierarchy");

    declareSimulator<ISimulator>(
        py::class_<ISimulator, std::shared_ptr<ISimulator>>(m, "ISimulator")
            .def("interp_order", &ISimulator::interporder));

    m.def("make_hierarchy", []() { return PHARE::amr::Hierarchy::make(); });
    m.def("make_simulator", [](std::shared_ptr<PHARE::amr::Hierarchy>& hier) {
        return std::shared_ptr<ISimulator>{std::move(PHARE::getSimulator(hier))};
    });

    py::class_<RuntimeDiagnosticInterface, std::shared_ptr<RuntimeDiagnosticInterface>>(
        m, "IDiagnosticsManager")
        .def("dump", &RuntimeDiagnosticInterface::dump, py::arg("timestamp"), py::arg("timestep"));

    m.def("make_diagnostic_manager", [](std::shared_ptr<ISimulator> const& sim,
                                        std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
        return std::make_shared<RuntimeDiagnosticInterface>(*sim, *hier);
    });

    m.def("mpi_size", []() {
        int mpi_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        return mpi_size;
    });
    m.def("mpi_rank", []() {
        int mpi_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        return mpi_rank;
    });

    declareDim<1>(m);
    declareDim<2>(m);
    declareDim<3>(m);

    core::apply(core::possibleSimulators(), [&](auto const& simType) { declare(m, simType); });

    declarePatchData<std::vector<double>, 1>(m, "PatchDataVectorDouble_1D");
    declarePatchData<std::vector<double>, 2>(m, "PatchDataVectorDouble_2D");
    declarePatchData<std::vector<double>, 3>(m, "PatchDataVectorDouble_3D");
}

} // namespace PHARE::pydata
