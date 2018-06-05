#ifndef PHARE_ION_POPULATION_H
#define PHARE_ION_POPULATION_H

#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>


#include "hybrid/hybrid_quantities.h"


namespace PHARE
{
template<typename ParticleArray, typename Field, typename VecField>
class IonPopulation
{
public:
    IonPopulation(std::string name, double mass)
        : name_{std::move(name)}
        , mass_{mass}
    {
    }

    double mass() const { return mass_; }

    std::string const& name() const { return name_; }




    bool isUsable() const
    {
        return domainParticles_ != nullptr && ghostParticles_ != nullptr
               && coarseToFineParticles_ != nullptr && rho_ != nullptr;
    }


    bool isSettable() const
    {
        return domainParticles_ == nullptr && ghostParticles_ == nullptr
               && coarseToFineParticles_ == nullptr && rho_ == nullptr;
    }



    ParticleArray& domainParticles()
    {
        if (isUsable())
        {
            return *domainParticles_;
        }
        else
        {
            throw std::runtime_error("Error - cannot provide access to particle buffers");
        }
    }



    ParticleArray& ghostParticles()
    {
        if (isUsable())
        {
            return *ghostParticles_;
        }
        else
        {
            throw std::runtime_error("Error - cannot provide access to particle buffers");
        }
    }



    ParticleArray& coarseToFineParticles()
    {
        if (isUsable())
        {
            return *coarseToFineParticles_;
        }
        else
        {
            throw std::runtime_error("Error - cannot provide access to particle buffers");
        }
    }


    Field const& density() const
    {
        if (isUsable())
        {
            return *rho_;
        }
        else
        {
            throw std::runtime_error("Error - cannot provide access to density field");
        }
    }


    Field& density()
    {
        return const_cast<Field&>(static_cast<const IonPopulation*>(this)->density);
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


    using field_type = Field;


    MomentProperties getFieldNamesAndQuantities() const
    {
        return {{{name_ + "_rho", HybridQuantity::Scalar::rho}}};
    }



    struct ParticleProperty
    {
        std::string name;
    };

    using ParticleProperties = std::vector<ParticleProperty>;


    std::vector<std::string> getParticleArrayNames() const
    {
        return {{{name_ + "_domain"}, {name_ + "_ghost"}, {name_ + "_coarseToFine"}}};
    }


    auto getSubResourcesObject() const { return std::forward_as_tuple(flux); }


    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



private:
    std::string name_;
    double mass_;

    Field* rho_{nullptr};
    VecField flux;

    ParticleArray* domainParticles_{nullptr};
    ParticleArray* ghostParticles_{nullptr};
    ParticleArray* coarseToFineParticles_{nullptr};
};

} // namespace PHARE

#endif
