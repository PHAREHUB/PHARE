#ifndef PHARE_TEST_CORE_DATA_MHDSTATE_MHDSTATE_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_MHDSTATE_MHDSTATE_FIXTURES_HPP

#include "core/models/mhd_state.hpp"
#include "core/models/quantities/mhd_quantities.hpp"

#include "initializer/data_provider.hpp"

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
        auto&& [_rho, _V, _B1, _P, _rhoV, _Etot1, _J, _E]
            = Super::getCompileTimeResourcesViewList();
        _rho.setBuffer(&rho);
        V.set_on(_V);
        B1.set_on(_B1);
        _P.setBuffer(&P);
        rhoV.set_on(_rhoV);
        _Etot1.setBuffer(&Etot1);
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
        , B1{dict["name"].template to<std::string>() + "_" + "B1", layout, MHDQuantity::Vector::B1}
        , P{dict["name"].template to<std::string>() + "_" + "P", layout, MHDQuantity::Scalar::P}
        , rhoV{dict["name"].template to<std::string>() + "_" + "rhoV", layout,
               MHDQuantity::Vector::rhoV}
        , Etot1{dict["name"].template to<std::string>() + "_" + "Etot1", layout,
                MHDQuantity::Scalar::Etot1}
        , J{dict["name"].template to<std::string>() + "_" + "J", layout, MHDQuantity::Vector::J}
        , E{dict["name"].template to<std::string>() + "_" + "E", layout, MHDQuantity::Vector::E}
    {
        _set();
    }

    template<typename GridLayout>
    UsableMHDState(GridLayout const& layout, std::string name)
        : Super{name}
        , rho{name + "_rho", layout, MHDQuantity::Scalar::rho}
        , V{name + "_V", layout, MHDQuantity::Vector::V}
        , B1{name + "_B1", layout, MHDQuantity::Vector::B1}
        , P{name + "_P", layout, MHDQuantity::Scalar::P}
        , rhoV{name + "_rhoV", layout, MHDQuantity::Vector::rhoV}
        , Etot1{name + "_Etot1", layout, MHDQuantity::Scalar::Etot1}
        , J{name + "_J", layout, MHDQuantity::Vector::J}
        , E{name + "_E", layout, MHDQuantity::Vector::E}
    {
        _set();
    }

    UsableMHDState(UsableMHDState const&) = delete;

    UsableMHDState(UsableMHDState&& that)
        : Super{std::forward<Super>(that)}
        , rho{std::move(that.rho)}
        , V{std::move(that.V)}
        , B1{std::move(that.B1)}
        , P{std::move(that.P)}
        , rhoV{std::move(that.rhoV)}
        , Etot1{std::move(that.Etot1)}
        , J{std::move(that.J)}
        , E{std::move(that.E)}
    {
        _set();
    }

    void set_on(Super& state)
    {
        auto&& [_rho, _V, _B1, _P, _rhoV, _Etot1, _J, _E] = state.getCompileTimeResourcesViewList();

        _rho.setBuffer(&rho);
        V.set_on(_V);
        B1.set_on(_B1);
        _P.setBuffer(&P);
        rhoV.set_on(_rhoV);
        _Etot1.setBuffer(&Etot1);
        J.set_on(_J);
        E.set_on(_E);
    }

    Super& super() { return *this; }
    Super const& super() const { return *this; }
    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }

    Grid_t rho;
    UsableVecFieldMHD<dim> V;
    UsableVecFieldMHD<dim> B1;
    Grid_t P;

    UsableVecFieldMHD<dim> rhoV;
    Grid_t Etot1;

    UsableVecFieldMHD<dim> J, E;
};

} // namespace PHARE::core

#endif
