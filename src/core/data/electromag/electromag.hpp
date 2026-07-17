#ifndef PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_HPP
#define PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_HPP


#include "core/data/grid/grid_tiles.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/def.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"

#include "core/utilities/types.hpp"
#include "initializer/data_provider.hpp"


#include <tuple>
#include <string>

namespace PHARE::core::basic
{


template<typename VecFieldT>
class Electromag
{
public:
    using vecfield_type                    = VecFieldT;
    static constexpr std::size_t dimension = VecFieldT::dimension;

    template<typename V>
    auto as(auto&& a, auto&&... args) _PHARE_ALL_FN_
    {
        return V{{a(E, args...)}, {a(B, args...)}};
    }

    template<typename V>
    auto as(auto&& a, auto&&... args) const _PHARE_ALL_FN_
    {
        return V{{a(E, args...)}, {a(B, args...)}};
    }

    VecFieldT E, B;
};


} // namespace PHARE::core::basic


namespace PHARE
{
namespace core
{
    template<typename VecFieldT>
    class Electromag : public basic::Electromag<VecFieldT>
    {
    public:
        using Super = basic::Electromag<VecFieldT>;
        using Super::B;
        using Super::E;
        using vecfield_type                    = VecFieldT;
        static constexpr std::size_t dimension = VecFieldT::dimension;

        explicit Electromag(std::string name)
            : Super{{name + "_E", HybridQuantity::Vector::E},
                    {name + "_B", HybridQuantity::Vector::B}}
            , Binit_{}
        {
        }

        explicit Electromag(initializer::PHAREDict const& dict)
            : Super{{dict["name"].template to<std::string>() + "_"
                         + dict["electric"]["name"].template to<std::string>(),
                     HybridQuantity::Vector::E},
                    {dict["name"].template to<std::string>() + "_"
                         + dict["magnetic"]["name"].template to<std::string>(),
                     HybridQuantity::Vector::B}}
            , Binit_{dict["magnetic"]["initializer"]}
        {
        }


        template<typename GridLayout>
        void initialize(GridLayout const& layout)
        {
            assert(B.isUsable());
            Binit_.initialize(B, layout);
        }


        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------

        NO_DISCARD bool isUsable() const { return E.isUsable() && B.isUsable(); }

        NO_DISCARD bool isSettable() const { return E.isSettable() && B.isSettable(); }

        NO_DISCARD auto getCompileTimeResourcesViewList() const
        {
            return std::forward_as_tuple(E, B);
        }

        NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(E, B); }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        void copyData(Electromag const& source)
        {
            E.copyData(source.E);
            B.copyData(source.B);
        }


        auto operator()() { return std::forward_as_tuple(E, B); }
        auto operator()() const { return std::forward_as_tuple(E, B); }




        Super& super() { return *this; }
        Super const& super() const { return *this; }
        auto& operator*() { return super(); }
        auto& operator*() const { return super(); }


    private:
        VecFieldInitializer<dimension> Binit_;
    };


} // namespace core
} // namespace PHARE


namespace PHARE::core
{

void check_electromag(auto const& em) _PHARE_ALL_FN_
{
#if PHARE_DEBUG && !PHARE_HAVE_GPU
    check_tensor_field(em.E);
    check_tensor_field(em.B);
#endif
}

} // namespace PHARE::core


#endif
