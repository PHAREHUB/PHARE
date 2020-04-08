
#include "phare/include.h"

#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;

namespace PHARE
{
template<size_t dim_, size_t interp_>
struct ContiguousParticles : public core::ContiguousParticles<dim_>
{
    static const int32_t babies = 2;
    static const size_t dim     = dim_;
    static const size_t interp  = interp_;
    using Super                 = core::ContiguousParticles<dim>;
    using Super::charge;
    using Super::delta;
    using Super::iCell;
    using Super::size;
    using Super::v;
    using Super::weight;

    ContiguousParticles(size_t s)
        : core::ContiguousParticles<dim>{s}
        , refineFactor{arr()}
        , splitter{refineFactor, babies}
    {
    }

    static constexpr std::array<int32_t, dim> arr()
    {
        if constexpr (dim == 1)
            return std::array<int32_t, dim>{babies};
        if constexpr (dim == 2)
            return std::array<int32_t, dim>{babies, babies};
        if constexpr (dim == 3)
            return std::array<int32_t, dim>{babies, babies, babies};
    }

    core::Particle<dim> particle(size_t i) const
    {
        core::Particle<dim> p;
        auto copy = [i](auto& from, auto& to, auto len) {
            using T = typename std::decay_t<decltype(from)>::value_type;
            memcpy(to.data(), &from[i * len], sizeof(T) * len);
        };
        copy(iCell, p.iCell, dim);
        copy(delta, p.delta, dim);
        copy(v, p.v, 3);
        p.weight = weight[i];
        p.charge = charge[i];
        return p;
    }

    void add(core::Particle<dim>& p, size_t i)
    {
        auto copy = [i](auto& from, auto& to, auto len) {
            using T = typename std::decay_t<decltype(from)>::value_type;
            memcpy(&to[i * len], from.data(), sizeof(T) * len);
        };
        copy(p.iCell, iCell, dim);
        copy(p.delta, delta, dim);
        copy(p.v, v, 3);
        weight[i] = p.weight;
        charge[i] = p.charge;
    }

    ContiguousParticles<dim, interp> split() const
    {
        ContiguousParticles<dim, interp> kinder(size() * babies);
        for (size_t i = 0; i < size(); i++)
        {
            std::vector<core::Particle<dim>> refinedParticles;
            auto part = particle(i);
            auto fine = amr::CoarseParticlePrefiner<dim>::prefine(part, refineFactor);
            splitter(fine, refinedParticles);
            for (size_t j = 0; j < babies; j++)
                kinder.add(refinedParticles[j], (i * babies) + j);
        }
        return kinder;
    }

    std::array<int32_t, dim> refineFactor;
    PHARE::amr::Split<dim, interp> splitter;
};


template<size_t dim, size_t interp>
void declareParticlesDim(py::module& m)
{
    using T          = ContiguousParticles<dim, interp>;
    std::string name = "ContiguousParticles_" + std::to_string(dim) + "_" + std::to_string(interp);
    py::class_<T, std::shared_ptr<T>>(m, name.c_str())
        .def(py::init<size_t>())
        .def_readwrite("iCell", &T::iCell)
        .def_readwrite("delta", &T::delta)
        .def_readwrite("weight", &T::weight)
        .def_readwrite("charge", &T::charge)
        .def_readwrite("v", &T::v)
        .def("size", &T::size)
        .def("split", &T::split);
}

template<size_t interp>
void declareParticles(py::module& m)
{
    declareParticlesDim<1, interp>(m);
    // declareParticlesDim<2, interp>(m);
    // declareParticlesDim<3, interp>(m);
}

PYBIND11_MODULE(test_simulator, m)
{
    declareParticles<1>(m);
    declareParticles<2>(m);
    declareParticles<3>(m);
}
} // namespace PHARE
