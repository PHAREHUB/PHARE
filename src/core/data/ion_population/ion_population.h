#ifndef PHARE_ION_POPULATION_H
#define PHARE_ION_POPULATION_H

#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>


#include "hybrid/hybrid_quantities.h"
#include "particle_pack.h"


namespace PHARE
{
template<typename ParticleArray, typename VecField>
class IonPopulation
{
public:
    IonPopulation(std::string name, double mass)
        : name_{std::move(name)}
        , mass_{mass}
        , flux_{"flux", HybridQuantity::Vector::V}
    {
    }

    using field_type                       = typename VecField::field_type;
    static constexpr std::size_t dimension = VecField::dimension;
    using particle_array_type              = ParticleArray;
    using particle_resource_type           = ParticlesPack<ParticleArray>;

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


    field_type const& density() const
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


    field_type& density()
    {
        return const_cast<field_type&>(static_cast<const IonPopulation*>(this)->density);
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



    struct ParticleProperty
    {
        std::string name;
    };

    using ParticleProperties = std::vector<ParticleProperty>;


    ParticleProperties getParticleArrayNames() const
    {
        return {{{name_ + "_domain"}, {name_ + "_ghost"}, {name_ + "_coarseToFine"}}};
    }


    void setBuffer(std::string const& bufferName, ParticlesPack<ParticleArray>* pack)
    {
        if (pack != nullptr)
        {
            domainParticles_       = pack->domainParticles;
            ghostParticles_        = pack->ghostParticles;
            coarseToFineParticles_ = pack->coarseToFineParticles;
        }
        else
        {
            domainParticles_       = nullptr;
            ghostParticles_        = nullptr;
            coarseToFineParticles_ = nullptr;
        }
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



    auto getSubResourcesObject() { return std::forward_as_tuple(flux_); }


    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



private:
    std::string name_;
    double mass_;
    VecField flux_;

    field_type* rho_{nullptr};

    ParticleArray* domainParticles_{nullptr};
    ParticleArray* ghostParticles_{nullptr};
    ParticleArray* coarseToFineParticles_{nullptr};
};

} // namespace PHARE

#endif
