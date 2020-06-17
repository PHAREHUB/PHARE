#include "phare/include.h"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "phare/include.h"
#include "core/utilities/types.h"
#include "core/utilities/mpi_utils.h"
#include "core/data/particles/particle_packer.h"




namespace py = pybind11;



namespace PHARE::pydata
{
template<typename Data>
struct PatchData
{
    std::string patchID;
    std::string origin;
    std::string cells;
    size_t nGhosts;
    Data data;

    PatchData() {}

    template<typename... Args>
    PatchData(Args&&... args)
        : data{std::forward<Args...>(args...)}
    {
    }
};


template<typename PatchData>
void setPatchData(PatchData& data, std::string patchID, std::string origin, std::string cells)
{
    data.patchID = patchID;
    data.origin  = origin;
    data.cells   = cells;
}

template<typename PatchData, typename GridLayout>
void setPatchDataFromGrid(PatchData& data, GridLayout& grid, std::string patchID)
{
    setPatchData(data, patchID, grid.origin().str(),
                 core::Point<uint32_t, GridLayout::dimension>{grid.nbrCells()}.str());
}

template<typename PatchData, typename Field, typename GridLayout>
void setPatchDataFromField(PatchData& data, Field const& field, GridLayout& grid,
                           std::string patchID)
{
    setPatchDataFromGrid(data, grid, patchID);
    data.nGhosts = static_cast<size_t>(
        GridLayout::nbrGhosts(GridLayout::centering(field.physicalQuantity())[0]));
    data.data.assign(field.data(), field.data() + field.size());
}


template<typename DataWrangler>
class PatchLevel
{
public:
    static constexpr size_t dimension     = DataWrangler::dimension;
    static constexpr size_t interp_order  = DataWrangler::interp_order;
    static constexpr size_t nbRefinedPart = DataWrangler::nbRefinedPart;

    using PHARETypes  = PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using HybridModel = typename PHARETypes::HybridModel_t;
    using GridLayout  = typename HybridModel::gridlayout_type;

    PatchLevel(amr::Hierarchy& hierarchy, HybridModel& model, size_t lvl)
        : lvl_(lvl)
        , hierarchy_{hierarchy}
        , model_{model}
    {
    }

    auto getDensity()
    {
        std::vector<PatchData<std::vector<double>>> patch_data;
        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, size_t /*iLevel*/) {
            setPatchDataFromField(patch_data.emplace_back(), ions.density(), grid, patchID);
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, ions);

        return patch_data;
    }

    auto getPopDensities()
    {
        using Inner = decltype(getDensity());

        std::unordered_map<std::string, Inner> pop_data;
        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, size_t /*iLevel*/) {
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

        auto visit = [&](GridLayout& grid, std::string patchID, size_t /*iLevel*/) {
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

    auto getPopFluxs()
    {
        using Inner = decltype(getPopDensities());

        std::unordered_map<std::string, Inner> pop_data;

        auto& ions = model_.state.ions;

        auto visit = [&](GridLayout& grid, std::string patchID, size_t /*iLevel*/) {
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

        auto visit = [&](GridLayout& grid, std::string patchID, size_t /*iLevel*/) {
            for (auto& vecFieldPtr : {&em.B, &em.E})
                getVecFields(*vecFieldPtr, em_data, grid, patchID, vecFieldPtr->name());
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, em);

        return em_data;
    }

    auto getParticles()
    {
        using Nested = std::vector<PatchData<core::ContiguousParticles<dimension>>>;
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

        auto visit = [&](GridLayout& grid, std::string patchID, size_t /*iLevel*/) {
            for (auto& pop : ions)
            {
                if (!pop_particles.count(pop.name()))
                    pop_particles.emplace(pop.name(), Inner());

                auto& inner = pop_particles.at(pop.name());

                getParticleData(inner, grid, patchID, "domain", pop.domainParticles());
                getParticleData(inner, grid, patchID, "patchGhost", pop.patchGhostParticles());
                getParticleData(inner, grid, patchID, "levelGhost", pop.levelGhostParticles());
            }
        };

        PHARE::amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(lvl_),
                                           *model_.resourcesManager, visit, ions);

        return pop_particles;
    }

private:
    amr::Hierarchy& hierarchy_;
    HybridModel& model_;
    size_t lvl_;
};

template<std::size_t _dimension, std::size_t _interp_order, size_t _nbRefinedPart>
class DataWrangler
{
public:
    using This                            = DataWrangler;
    static constexpr size_t dimension     = _dimension;
    static constexpr size_t interp_order  = _interp_order;
    static constexpr size_t nbRefinedPart = _nbRefinedPart;

    using PHARETypes  = PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using HybridModel = typename PHARETypes::HybridModel_t;

    DataWrangler(std::shared_ptr<ISimulator> const& simulator,
                 std::shared_ptr<amr::Hierarchy> const& hierarchy)
        : simulator_{simulator}
        , hierarchy_{hierarchy}
    {
        auto dict          = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
        auto dim           = dict["simulation"]["dimension"].template to<int>();
        auto interpOrder   = dict["simulation"]["interp_order"].template to<int>();
        auto nbRefinedPart = dict["simulation"]["refined_particle_nbr"].template to<int>();
        if (!core::makeAtRuntime<Maker>(dim, interpOrder, nbRefinedPart, Maker{*this}))
            throw std::runtime_error("Runtime diagnostic deduction failed");
    }

    auto getPatchLevel(size_t lvl)
    {
        return PatchLevel<This>{*hierarchy_, *simulator_ptr_->getHybridModel(), lvl};
    }

    auto sort_merge_1d(std::vector<PatchData<std::vector<double>>> const&& input,
                       bool shared_patch_border = false)
    {
        std::vector<std::pair<double, const PatchData<std::vector<double>>*>> sorted;
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

    auto sync(std::vector<PatchData<std::vector<double>>> const& input)
    {
        int mpi_size = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        std::vector<PatchData<std::vector<double>>> collected;

        auto collect = [&](auto& patch_data) {
            auto patchIDs = core::mpi::collect(patch_data.patchID, mpi_size);
            auto origins  = core::mpi::collect(patch_data.origin, mpi_size);
            auto pCells   = core::mpi::collect(patch_data.cells, mpi_size);
            auto ghosts   = core::mpi::collect(patch_data.nGhosts, mpi_size);
            auto datas    = core::mpi::collect(patch_data.data, mpi_size);

            for (int i = 0; i < mpi_size; i++)
            {
                auto& data = collected.emplace_back();
                setPatchData(data, patchIDs[i], origins[i], pCells[i]);
                data.nGhosts = ghosts[i];
                data.data    = std::move(datas[i]);
            }
        };

        auto max = core::mpi::max(input.size(), mpi_size);

        PatchData<std::vector<double>> empty;

        for (size_t i = 0; i < max; i++)
        {
            if (i < input.size())
                collect(input[i]);
            else
                collect(empty);
        }
        return collected;
    }

    auto sync_merge(std::vector<PatchData<std::vector<double>>> const& input, bool primal)
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
        bool operator()(std::size_t userDim, std::size_t userInterpOrder, size_t userNbRefinedPart,
                        Dimension dimension_fn, InterpOrder interp_order_fn,
                        NbRefinedPart nbRefinedPart_fn)
        {
            if (userDim == dimension_fn() and userInterpOrder == interp_order_fn()
                and userNbRefinedPart == nbRefinedPart_fn())
            {
                size_t constexpr d  = dimension_fn();
                size_t constexpr io = interp_order_fn();
                size_t constexpr nb = nbRefinedPart_fn();

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



template<typename Type>
void declarePatchData(py::module& m, std::string key)
{
    using PatchDataType = PatchData<Type>;
    py::class_<PatchDataType>(m, key.c_str())
        .def_readonly("patchID", &PatchDataType::patchID)
        .def_readonly("origin", &PatchDataType::origin)
        .def_readonly("cells", &PatchDataType::cells)
        .def_readonly("nGhosts", &PatchDataType::nGhosts)
        .def_readonly("data", &PatchDataType::data);
}

template<size_t dim>
void declareDim(py::module& m)
{
    using CP         = core::ContiguousParticles<dim>;
    std::string name = "ContiguousParticles_" + std::to_string(dim);
    py::class_<CP, std::shared_ptr<CP>>(m, name.c_str())
        .def(py::init<size_t>())
        .def_readwrite("iCell", &CP::iCell)
        .def_readwrite("delta", &CP::delta)
        .def_readwrite("weight", &CP::weight)
        .def_readwrite("charge", &CP::charge)
        .def_readwrite("v", &CP::v)
        .def("size", &CP::size);

    name = "PatchData" + name;
    declarePatchData<CP>(m, name.c_str());
}



template<typename _dim, typename _interp, typename _nbRefinedPart>
void declare(py::module& m)
{
    constexpr auto dim           = _dim{}();
    constexpr auto interp        = _interp{}();
    constexpr auto nbRefinedPart = _nbRefinedPart{}();


    std::string type_string = "_" + std::to_string(dim) + "_" + std::to_string(interp) + "_"
                              + std::to_string(nbRefinedPart);

    using DW         = DataWrangler<dim, interp, nbRefinedPart>;
    std::string name = "DataWrangler" + type_string;
    py::class_<DW, std::shared_ptr<DW>>(m, name.c_str())
        .def(py::init<std::shared_ptr<ISimulator> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def("sync_merge", &DW::sync_merge)
        .def("getPatchLevel", &DW::getPatchLevel);

    using PL = PatchLevel<DW>;
    name     = "PatchLevel" + type_string;
    py::class_<PL, std::shared_ptr<PL>>(m, name.c_str())
        .def("getEM", &PL::getEM)
        .def("getDensity", &PL::getDensity)
        .def("getBulkVelocity", &PL::getBulkVelocity)
        .def("getPopDensities", &PL::getPopDensities)
        .def("getPopFluxs", &PL::getPopFluxs)
        .def("getParticles", &PL::getParticles);

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

class StaticSamraiLifeCycle : public SamraiLifeCycle
{
public:
    inline static StaticSamraiLifeCycle& INSTANCE()
    {
        static StaticSamraiLifeCycle i;
        return i;
    }
};

PYBIND11_MODULE(cpp, m)
{
    StaticSamraiLifeCycle::INSTANCE(); // init

    py::class_<PHARE::amr::Hierarchy, std::shared_ptr<PHARE::amr::Hierarchy>>(m, "AMRHierarchy");

    py::class_<ISimulator, std::shared_ptr<ISimulator>>(m, "ISimulator")
        .def("initialize", &PHARE::ISimulator::initialize)
        .def("advance", &PHARE::ISimulator::advance)
        .def("startTime", &PHARE::ISimulator::startTime)
        .def("currentTime", &PHARE::ISimulator::currentTime)
        .def("endTime", &PHARE::ISimulator::endTime)
        .def("timeStep", &PHARE::ISimulator::timeStep)
        .def("to_str", &PHARE::ISimulator::to_str);

    m.def("make_hierarchy", []() { return PHARE::amr::Hierarchy::make(); });
    m.def("make_simulator", [](std::shared_ptr<PHARE::amr::Hierarchy>& hier) {
        return std::shared_ptr<ISimulator>{std::move(PHARE::getSimulator(hier))};
    });

    py::class_<RuntimeDiagnosticInterface, std::shared_ptr<RuntimeDiagnosticInterface>>(
        m, "IDiagnosticsManager")
        .def("dump", &RuntimeDiagnosticInterface::dump);

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
    m.def("reset", []() {
        py::gil_scoped_release release;
        StaticSamraiLifeCycle::reset();
    });

    declareDim<1>(m);
    declareDim<2>(m);
    declareDim<3>(m);

    core::apply(core::possibleSimulators(), [&](auto const& simType) { declare(m, simType); });

    declarePatchData<std::vector<double>>(m, "PatchDataVectorDouble");
}

} // namespace PHARE::pydata
