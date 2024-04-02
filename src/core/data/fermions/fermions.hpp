#ifndef PHARE_FERMIONS_HPP
#define PHARE_FERMIONS_HPP

#include "core/data/ions/ions.hpp"
#include "core/data/fermions/pic_electrons.hpp"

namespace PHARE::core
{
    template<typename Ions, typename PICElectrons>
    class Fermions
    {
    public:
        static constexpr auto dimension = Ions::dimension;

        explicit Fermions(PHARE::initializer::PHAREDict const& dict)
        : ions{dict["ions"]}
        , electrons{} // TODO

        Ions ions;
        PICElectrons electrons;

        // ResourcesUser interface??
        
    };
} // namespace PHARE::core

#endif
