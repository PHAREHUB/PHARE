#ifndef PHARE_TEST_CORE_DATA_MHDSTATE_MHDSTATE_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_MHDSTATE_MHDSTATE_FIXTURES_HPP

#include "core/mhd/mhd_quantities.hpp"
#include "core/models/mhd_state.hpp"
#include "initializer/data_provider.hpp"
#include "phare_core.hpp"
#include "tests/core/data/field/test_field_fixtures_mhd.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

namespace PHARE::core
{
using namespace PHARE::initializer;

template<std::size_t dim>
class UsableMHDState : public MHDState<VecFieldMHD<dim>>
{
    using Super = MHDState<VecFieldMHD<dim>>;

    void _set()
    {
        auto&& [_rho, _V, _B, _P, _rhoV, _Etot, _J, _E] = Super::getCompileTimeResourcesViewList();
        _rho.setBuffer(&rho);
        V.set_on(_V);
        B.set_on(_B);
        _P.setBuffer(&P);
        rhoV.set_on(_rhoV);
        _Etot.setBuffer(&Etot);
        J.set_on(_J);
        E.set_on(_E);
    }

public:
    using Array_t = NdArrayVector<dim, double, /*c_ordering*/ true>;
    using Grid_t  = Grid<Array_t, MHDQuantity::Scalar>;

    template<typename GridLayout>
    UsableMHDState(GridLayout const& layout, PHAREDict const& dict)
        : Super{dict}
        , rho{dict["name"].template to<std::string>() + "_" + "rho", layout,
              MHDQuantity::Scalar::rho}
        , V{dict["name"].template to<std::string>() + "_" + "V", layout, MHDQuantity::Vector::V}
        , B{dict["name"].template to<std::string>() + "_" + "B", layout, MHDQuantity::Vector::B}
        , P{dict["name"].template to<std::string>() + "_" + "P", layout, MHDQuantity::Scalar::P}
        , rhoV{dict["name"].template to<std::string>() + "_" + "rhoV", layout,
               MHDQuantity::Vector::rhoV}
        , Etot{dict["name"].template to<std::string>() + "_" + "Etot", layout,
               MHDQuantity::Scalar::Etot}
        , J{dict["name"].template to<std::string>() + "_" + "J", layout, MHDQuantity::Vector::J}
        , E{dict["name"].template to<std::string>() + "_" + "E", layout, MHDQuantity::Vector::E}
    {
        _set();
    }

    UsableMHDState(UsableMHDState const&) = delete;

    UsableMHDState(UsableMHDState&& that)
        : Super{std::forward<Super>(that)}
        , rho{std::move(that.rho)}
        , V{std::move(that.V)}
        , B{std::move(that.B)}
        , P{std::move(that.P)}
        , rhoV{std::move(that.rhoV)}
        , Etot{std::move(that.Etot)}
        , J{std::move(that.J)}
        , E{std::move(that.E)}
    {
        _set();
    }

    Super& super() { return *this; }
    Super const& super() const { return *this; }
    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }

    Grid_t rho;
    UsableVecFieldMHD<dim> V;
    UsableVecFieldMHD<dim> B;
    Grid_t P;

    UsableVecFieldMHD<dim> rhoV;
    Grid_t Etot;

    UsableVecFieldMHD<dim> J, E;
};

} // namespace PHARE::core

#endif
