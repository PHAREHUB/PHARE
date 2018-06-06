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
        , bulkVelocity_{name_ + "_bulkVel", HybridQuantity::Vector::V}
    {
    }


    auto begin() { return std::begin(populations_); }
    auto end() { return std::end(populations_); }


    bool isUsable() const
    {
        bool usable = rho_ != nullptr && bulkVelocity_.isUsable();
        for (auto const& pop : populations_)
        {
            usable = usable && pop.isUsable();
        }
        return usable;
    }



    bool isSettable() const
    {
        bool settable = rho_ == nullptr && bulkVelocity_.isSettable();
        for (auto const& pop : populations_)
        {
            settable = settable && pop.isSettable();
        }
        return settable;
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------


    struct MomentsProperty
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using MomentProperties = std::vector<MomentsProperty>;

    MomentProperties getFieldNamesAndQuantities() const
    {
        return {{{name_ + "_rho", HybridQuantity::Scalar::rho}}};
    }



    void setBuffer(std::string const& bufferName, field_type* field)
    {
        if (bufferName == name_ + "_rho")
        {
            rho_ = field;
        }
        else
        {
            throw std::runtime_error("Error - invalid density buffer name");
        }
    }



    std::vector<IonPopulation>& getRunTimeResourcesUserList() { return populations_; }

    auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(bulkVelocity_); }



    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

private:
    std::string name_;
    field_type* rho_{nullptr};
    vecfield_type bulkVelocity_;
    std::vector<IonPopulation> populations_; // TODO we have to name this so they are unique
                                             // although only 1 Ions should exist.
};
} // namespace PHARE

#endif
