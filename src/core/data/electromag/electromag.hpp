#ifndef PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_HPP
#define PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_HPP

#include <string>
#include <tuple>

#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"
#include "initializer/data_provider.hpp"
#include "core/def.hpp"


namespace PHARE
{
namespace core
{
    template<typename VecFieldT>
    class Electromag
    {
    public:
        static constexpr std::size_t dimension = VecFieldT::dimension;

        explicit Electromag(std::string name)
            : E{name + "_E", HybridQuantity::Vector::E}
            , B{name + "_B", HybridQuantity::Vector::B}
            , Binit_{}
        {
        }

        explicit Electromag(initializer::PHAREDict const& dict)
            : E{dict["name"].template to<std::string>() + "_"
                    + dict["electric"]["name"].template to<std::string>(),
                HybridQuantity::Vector::E}
            , B{dict["name"].template to<std::string>() + "_"
                    + dict["magnetic"]["name"].template to<std::string>(),
                HybridQuantity::Vector::B}
            , Binit_{dict["magnetic"]["initializer"]}
        {
        }

        using vecfield_type = VecFieldT;


        template<typename GridLayout>
        void initialize(GridLayout const& layout)
        {
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

        VecFieldT E;
        VecFieldT B;

    private:
        VecFieldInitializer<dimension> Binit_;
    };
} // namespace core
} // namespace PHARE
#endif
