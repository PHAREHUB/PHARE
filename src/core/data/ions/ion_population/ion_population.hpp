#ifndef PHARE_ION_POPULATION_HPP
#define PHARE_ION_POPULATION_HPP

#include <tuple>
#include <algorithm>
#include <string>
#include <sstream>


#include "core/def.hpp"
#include "particle_pack.hpp"
#include "initializer/data_provider.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

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
            , charge_{initializer["particle_initializer"]["charge"].template to<double>()}
            , flux_{name_ + "_flux", HybridQuantity::Vector::V}
            , momentumTensor_{name_ + "_momentumTensor", HybridQuantity::Tensor::M}
            , particleDensity_{name_ + "_particleDensity", HybridQuantity::Scalar::rho}
            , chargeDensity_{name_ + "_chargeDensity", HybridQuantity::Scalar::rho}
            , particles_{name_}
            , particleInitializerInfo_{initializer["particle_initializer"]}
        {
        }

        IonPopulation(IonPopulation const&)            = default;
        IonPopulation(IonPopulation&&)                 = default;
        IonPopulation& operator=(IonPopulation const&) = default;
        IonPopulation& operator=(IonPopulation&&)      = default;

        NO_DISCARD double mass() const { return mass_; }

        //! all particles of a population carry this charge, see particlesHaveCharge()
        NO_DISCARD double charge() const { return charge_; }

        NO_DISCARD std::string const& name() const { return name_; }



        NO_DISCARD auto const& particleInitializerInfo() const
        {
            assert(particleInitializerInfo_.contains("density"));
            return particleInitializerInfo_;
        }


        NO_DISCARD bool isUsable() const
        {
            return core::isUsable(particles_, particleDensity_, chargeDensity_, flux_,
                                  momentumTensor_);
        }

        NO_DISCARD bool isSettable() const
        {
            return core::isSettable(particles_, particleDensity_, chargeDensity_, flux_,
                                    momentumTensor_);
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

        NO_DISCARD field_type const& particleDensity() const { return particleDensity_; }
        NO_DISCARD field_type& particleDensity() { return particleDensity_; }

        NO_DISCARD field_type const& chargeDensity() const { return chargeDensity_; }
        NO_DISCARD field_type& chargeDensity() { return chargeDensity_; }

        NO_DISCARD VecField const& flux() const { return flux_; }
        NO_DISCARD VecField& flux() { return flux_; }

        NO_DISCARD TensorField const& momentumTensor() const { return momentumTensor_; }
        NO_DISCARD TensorField& momentumTensor() { return momentumTensor_; }


        //! true if every particle, in every array set on the population, carries charge()
        NO_DISCARD bool particlesHaveCharge() const
        {
            auto const haveCharge = [this](ParticleArray const* particles) {
                return particles == nullptr
                       or std::all_of(std::begin(*particles), std::end(*particles),
                                      [this](auto const& part) { return part.charge == charge_; });
            };
            return haveCharge(particles_._domainParticles)
                   and haveCharge(particles_._patchGhostParticles)
                   and haveCharge(particles_._levelGhostParticles)
                   and haveCharge(particles_._levelGhostParticlesOld)
                   and haveCharge(particles_._levelGhostParticlesNew);
        }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------




        NO_DISCARD auto getCompileTimeResourcesViewList()
        {
            return std::forward_as_tuple(flux_, momentumTensor_, particleDensity_, chargeDensity_,
                                         particles_);
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
        double charge_;
        VecField flux_;
        TensorField momentumTensor_;
        field_type particleDensity_;
        field_type chargeDensity_;
        ParticlesPack<ParticleArray> particles_;
        initializer::PHAREDict particleInitializerInfo_;
    };
} // namespace core
} // namespace PHARE

#endif
