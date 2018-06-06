#ifndef PHARE_IONS_H
#define PHARE_IONS_H

#include <iterator>

#include "data/ion_population/ion_population.h"
#include "hybrid/hybrid_quantities.h"

namespace PHARE
{
template<typename IonPopulation>
class Ions
{
public:
    using field_type    = typename IonPopulation::field_type;
    using vecfield_type = typename IonPopulation::vecfield_type;

    Ions(std::string name)
        : name_{std::move(name)}
    {
    }


    auto begin() { return std::begin(populations_); }
    auto end() { return std::end(populations_); }


    bool isUsable() const
    {
        bool usable = false;
        for (auto const& pop : populations_)
        {
            usable = usable && pop.isUsable();
        }
        return usable;
    }



    bool isSettable() const
    {
        bool settable = true;
        for (auto const& pop : populations_)
        {
            settable = settable && pop.isSettable();
        }
        return settable;
    }


    struct MomentsProperty
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using MomentProperties = std::vector<MomentsProperty>;

    MomentProperties getFieldnamesAndQuantities() const
    {
        return {{{name_ + "rho_", HybridQuantity::Scalar::rho}}};
    }



    auto getSubResourcesObjet()
    {
        //
    }



private:
    std::string name_;
    std::vector<IonPopulation> populations_;
};
} // namespace PHARE

#endif
