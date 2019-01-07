#ifndef PHARE_ION_POPULATION_H
#define PHARE_ION_POPULATION_H

#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>


#include "data/ions/particle_initializers/particle_initializer.h"
#include "hybrid/hybrid_quantities.h"
#include "particle_pack.h"


namespace PHARE
{
template<typename ParticleArray, typename VecField, typename GridLayout>
class IonPopulation
{
public:
    using field_type                       = typename VecField::field_type;
    static constexpr std::size_t dimension = VecField::dimension;
    using particle_array_type              = ParticleArray;
    using particle_resource_type           = ParticlesPack<ParticleArray>;
    using vecfield_type                    = VecField;
    using particle_initializer_type        = ParticleInitializer<ParticleArray, GridLayout>;


    IonPopulation(std::string name, double mass,
                  std::unique_ptr<particle_initializer_type> particleInitializer)
        : name_{std::move(name)}
        , mass_{mass}
        , flux_{name_ + "_flux", HybridQuantity::Vector::V}
        , particleInitializer_{std::move(particleInitializer)}
    {
    }

    double mass() const { return mass_; }

    std::string const& name() const { return name_; }




    bool isUsable() const { return particles_ != nullptr && rho_ != nullptr && flux_.isUsable(); }


    bool isSettable() const
    {
        return particles_ == nullptr && rho_ == nullptr && flux_.isSettable();
    }



    void loadParticles(GridLayout const& layout)
    {
        if (isUsable())
        {
            particleInitializer_->loadParticles(*particles_->domainParticles, layout);
        }
        else
        {
            throw std::runtime_error("Error - cannot load particles, IonPopulation not usable");
        }
    }




    auto nbrParticles() const
    {
        if (isUsable())
        {
            return particles_->domainParticles->size();
        }
        else
        {
            throw std::runtime_error("Error - cannot access to particles");
        }
    }




    ParticleArray& domainParticles()
    {
        if (isUsable())
        {
            return *particles_->domainParticles;
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
            return *particles_->ghostParticles;
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
            return *particles_->coarseToFineParticles;
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
        return const_cast<field_type&>(static_cast<const IonPopulation*>(this)->density());
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

    ParticleProperties getParticleArrayNames() const { return {{{name_}}}; }




    void setBuffer(std::string const& bufferName, ParticlesPack<ParticleArray>* pack)
    {
        if (bufferName == name_)
            particles_ = pack;
        else
            throw std::runtime_error("Error - invalid particle resource name");
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



    auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(flux_); }


    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



private:
    std::string name_;
    double mass_;
    VecField flux_;
    field_type* rho_{nullptr};
    ParticlesPack<ParticleArray>* particles_{nullptr};
    std::unique_ptr<particle_initializer_type> particleInitializer_;
};

} // namespace PHARE

#endif
