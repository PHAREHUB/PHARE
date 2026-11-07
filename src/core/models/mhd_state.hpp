#ifndef PHARE_MHD_STATE_HPP
#define PHARE_MHD_STATE_HPP

#include "core/def.hpp"
#include "core/models/physical_state.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"

#include "initializer/data_provider.hpp"

#include <string>
#include <sstream>

namespace PHARE::core::basic
{
// the fields of an MHDState without its initializers, a view (e.g. of one tile) is made via as
template<typename VecFieldT>
struct MHDState
{
    using vecfield_type = VecFieldT;
    using field_type    = VecFieldT::field_type;

    static constexpr auto dimension = VecFieldT::dimension;

    template<typename VF>
    using view_t = MHDState<VF>;

    template<typename View>
    auto as(auto&& a, auto&&... args)
    {
        return View{{a(rho, args...)},  {a(V, args...)},    {a(B, args...)}, {a(P, args...)},
                    {a(rhoV, args...)}, {a(Etot, args...)}, {a(E, args...)}, {a(J, args...)}};
    }

    template<typename View>
    auto as(auto&& a, auto&&... args) const
    {
        return View{{a(rho, args...)},  {a(V, args...)},    {a(B, args...)}, {a(P, args...)},
                    {a(rhoV, args...)}, {a(Etot, args...)}, {a(E, args...)}, {a(J, args...)}};
    }

    field_type rho;
    VecFieldT V;
    VecFieldT B;
    field_type P;

    VecFieldT rhoV;
    field_type Etot;

    VecFieldT E;
    VecFieldT J;
};

} // namespace PHARE::core::basic


namespace PHARE
{
namespace core
{

    template<typename VecFieldT>
    class MHDState : public basic::MHDState<VecFieldT>, public IPhysicalState
    {
    public:
        using Super         = basic::MHDState<VecFieldT>;
        using vecfield_type = VecFieldT;
        using field_type    = VecFieldT::field_type;

        using Super::B;
        using Super::E;
        using Super::Etot;
        using Super::J;
        using Super::P;
        using Super::rho;
        using Super::rhoV;
        using Super::V;

        static constexpr auto dimension = VecFieldT::dimension;

        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------

        NO_DISCARD bool isUsable() const
        {
            return rho.isUsable() and V.isUsable() and B.isUsable() and P.isUsable()
                   and rhoV.isUsable() and Etot.isUsable() and J.isUsable() and E.isUsable();
        }

        NO_DISCARD bool isSettable() const
        {
            return rho.isSettable() and V.isSettable() and B.isSettable() and P.isSettable()
                   and rhoV.isSettable() and Etot.isSettable() and J.isSettable()
                   and E.isSettable();
        }

        NO_DISCARD auto getCompileTimeResourcesViewList() const
        {
            return std::forward_as_tuple(rho, V, B, P, rhoV, Etot, J, E);
        }

        NO_DISCARD auto getCompileTimeResourcesViewList()
        {
            return std::forward_as_tuple(rho, V, B, P, rhoV, Etot, J, E);
        }

        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        NO_DISCARD std::string to_str() const
        {
            std::stringstream ss;
            ss << "MHD State\n";
            ss << "------------------------------------\n";
            return ss.str();
        }

        MHDState(PHARE::initializer::PHAREDict const& dict)
            : MHDState{dict["name"].template to<std::string>(), dict}
        {
        }

        MHDState(std::string const& name, PHARE::initializer::PHAREDict const& dict)
            : Super{{name + "_" + "rho", MHDQuantity::Scalar::rho},
                    {name + "_" + "V", MHDQuantity::Vector::V},
                    {name + "_" + "B", MHDQuantity::Vector::B},
                    {name + "_" + "P", MHDQuantity::Scalar::P},
                    {name + "_" + "rhoV", MHDQuantity::Vector::rhoV},
                    {name + "_" + "Etot", MHDQuantity::Scalar::Etot},
                    {name + "_" + "E", MHDQuantity::Vector::E},
                    {name + "_" + "J", MHDQuantity::Vector::J}}
            , rhoinit_{dict["density"]["initializer"]
                           .template to<initializer::InitFunction<dimension>>()}
            , Vinit_{dict["velocity"]["initializer"]}
            , Binit_{dict["magnetic"]["initializer"]}
            , Pinit_{dict["pressure"]["initializer"]
                         .template to<initializer::InitFunction<dimension>>()}
            , gamma_{dict["to_conservative_init"]["heat_capacity_ratio"].template to<double>()}
        {
        }

        MHDState(std::string name)
            : Super{{name + "_" + "rho", MHDQuantity::Scalar::rho},
                    {name + "_" + "V", MHDQuantity::Vector::V},
                    {name + "_" + "B", MHDQuantity::Vector::B},
                    {name + "_" + "P", MHDQuantity::Scalar::P},
                    {name + "_" + "rhoV", MHDQuantity::Vector::rhoV},
                    {name + "_" + "Etot", MHDQuantity::Scalar::Etot},
                    {name + "_" + "E", MHDQuantity::Vector::E},
                    {name + "_" + "J", MHDQuantity::Vector::J}}
            , gamma_{}
        {
        }

        template<typename GridLayout>
        void initialize(GridLayout const& layout)
        {
            FieldUserFunctionInitializer::initialize(rho, layout, rhoinit_);
            Vinit_.initialize(V, layout);
            Binit_.initialize(B, layout);
            FieldUserFunctionInitializer::initialize(P, layout, Pinit_);

            ToConservativeConverter{layout, gamma_}(
                rho, V, B, P, rhoV, Etot); // initial to conservative conversion because we
                                           // store conservative quantities on the grid
        }

    private:
        initializer::InitFunction<dimension> rhoinit_;
        VecFieldInitializer<dimension> Vinit_;
        VecFieldInitializer<dimension> Binit_;
        initializer::InitFunction<dimension> Pinit_;

        double const gamma_;
    };
} // namespace core
} // namespace PHARE

#endif // PHARE_MHD_STATE_HPP
