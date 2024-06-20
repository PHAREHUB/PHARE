#ifndef PHARE_ION_POPULATION_HPP
#define PHARE_ION_POPULATION_HPP

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <array>


#include "core/def.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "initializer/data_provider.hpp"
#include "particle_pack.hpp"

namespace PHARE
{
namespace core
{
    template<typename ParticleArray, typename VecField, typename TensorField>
    class IonPopulation
    {
    public:
        using field_type                       = typename VecField::field_type;
        static constexpr std::size_t dimension = VecField::dimension;
        using particle_array_type              = ParticleArray;
        using vecfield_type                    = VecField;
        using tensorfield_type                 = TensorField;


        IonPopulation(initializer::PHAREDict const& initializer)
            : name_{initializer["name"].template to<std::string>()}
            , mass_{initializer["mass"].template to<double>()}
            , flux_{name_ + "_flux", HybridQuantity::Vector::V}
            , momentumTensor_{name_ + "_momentumTensor", HybridQuantity::Tensor::M}
            , particleDensity_{name_ + "_particleDensity", HybridQuantity::Scalar::rho}
            // , chargeDensity_{name_ + "_chargeDensity", HybridQuantity::Scalar::rho}
            , particles_{name_}
            , particleInitializerInfo_{initializer["particle_initializer"]}
        {
        }


        NO_DISCARD double mass() const { return mass_; }

        NO_DISCARD std::string const& name() const { return name_; }

        NO_DISCARD auto const& particleInitializerInfo() const { return particleInitializerInfo_; }

        NO_DISCARD bool isUsable() const
        {
            // return particles_.isUsable() && particleDensity_.isUsable() && chargeDensity_.isUsable() && flux_.isUsable()
            return particles_.isUsable() && particleDensity_.isUsable() && flux_.isUsable()
                   && momentumTensor_.isUsable();
        }

        NO_DISCARD bool isSettable() const
        {
            // return particles_.isSettable() && particleDensity_.isSettable() && chargeDensity_.isSettable() && flux_.isSettable()
            return particles_.isSettable() && particleDensity_.isSettable() && flux_.isSettable()
                   && momentumTensor_.isSettable();
        }

        NO_DISCARD auto& domainParticles() const { return particles_.domainParticles(); }
        NO_DISCARD auto& domainParticles() { return particles_.domainParticles(); }

        NO_DISCARD auto& patchGhostParticles() const { return particles_.patchGhostParticles(); }
        NO_DISCARD auto& patchGhostParticles() { return particles_.patchGhostParticles(); }

        NO_DISCARD auto& levelGhostParticles() const { return particles_.levelGhostParticles(); }
        NO_DISCARD auto& levelGhostParticles() { return particles_.levelGhostParticles(); }

        NO_DISCARD auto& levelGhostParticlesOld() { return particles_.levelGhostParticlesOld(); }
        NO_DISCARD auto& levelGhostParticlesOld() const
        {
            return particles_.levelGhostParticlesOld();
        }

        NO_DISCARD auto& levelGhostParticlesNew() { return particles_.levelGhostParticlesNew(); }
        NO_DISCARD auto& levelGhostParticlesNew() const
        {
            return particles_.levelGhostParticlesNew();
        }

        NO_DISCARD field_type const& density() const { return particleDensity_; }
        NO_DISCARD field_type& density() { return particleDensity_; }

        NO_DISCARD field_type const& chargeDensity() const { return particleDensity_; } // TODO ouam : will return ctahgeDensity_
        NO_DISCARD field_type& chargeDensity() { return particleDensity_; } // TODO ouam : will return ctahgeDensity_

        NO_DISCARD VecField const& flux() const { return flux_; }
        NO_DISCARD VecField& flux() { return flux_; }

        NO_DISCARD TensorField const& momentumTensor() const { return momentumTensor_; }
        NO_DISCARD TensorField& momentumTensor() { return momentumTensor_; }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------




        NO_DISCARD auto getCompileTimeResourcesViewList()
        {
            // return std::forward_as_tuple(flux_, momentumTensor_, particleDensity_, chargeDensity_, particles_);
            return std::forward_as_tuple(flux_, momentumTensor_, particleDensity_, particles_);
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------



        NO_DISCARD std::string to_str()
        {
            std::stringstream ss;
            ss << "Ions Population\n";
            ss << "------------------------------------\n";
            ss << "name                : " << name() << "\n";
            return ss.str();
        }

    private:
        std::string name_;
        double mass_;
        VecField flux_;
        TensorField momentumTensor_;
        field_type particleDensity_;
        // field_type chargeDensity_;
        ParticlesPack<ParticleArray> particles_;
        initializer::PHAREDict const& particleInitializerInfo_;
    };
} // namespace core
} // namespace PHARE

#endif
