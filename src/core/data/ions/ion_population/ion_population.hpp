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
    template<typename ParticleArray, typename VecField, typename TensorField, typename GridLayout>
    class IonPopulation
    {
    public:
        using field_type                       = typename VecField::field_type;
        static constexpr std::size_t dimension = VecField::dimension;
        using particle_array_type              = ParticleArray;
        using particle_resource_type           = ParticlesPack<ParticleArray>;
        using vecfield_type                    = VecField;
        using tensorfield_type                 = TensorField;


        IonPopulation(initializer::PHAREDict const& initializer)
            : name_{initializer["name"].template to<std::string>()}
            , mass_{initializer["mass"].template to<double>()}
            , flux_{name_ + "_flux", HybridQuantity::Vector::V}
            , momentumTensor_{name_ + "_momentumTensor", HybridQuantity::Tensor::M}
            , particleInitializerInfo_{initializer["particle_initializer"]}
        {
        }


        NO_DISCARD double mass() const { return mass_; }

        NO_DISCARD std::string const& name() const { return name_; }


        NO_DISCARD auto const& particleInitializerInfo() const { return particleInitializerInfo_; }



        NO_DISCARD bool isUsable() const
        {
            return particles_ != nullptr && rho_ != nullptr && flux_.isUsable()
                   && momentumTensor_.isUsable();
        }


        NO_DISCARD bool isSettable() const
        {
            return particles_ == nullptr && rho_ == nullptr && flux_.isSettable()
                   && momentumTensor_.isSettable();
        }




        NO_DISCARD auto nbrParticles() const
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




        NO_DISCARD auto& domainParticles() const
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

        NO_DISCARD auto& domainParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->domainParticles());
        }



        NO_DISCARD auto& patchGhostParticles() const
        {
            if (isUsable())
            {
                return *particles_->patchGhostParticles;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }

        NO_DISCARD auto& patchGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->patchGhostParticles());
        }


        NO_DISCARD auto& levelGhostParticles() const
        {
            if (isUsable())
            {
                return *particles_->levelGhostParticles;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }

        NO_DISCARD auto& levelGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->levelGhostParticles());
        }




        NO_DISCARD ParticleArray& levelGhostParticlesOld()
        {
            if (isUsable())
            {
                return *particles_->levelGhostParticlesOld;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }



        NO_DISCARD ParticleArray& levelGhostParticlesNew()
        {
            if (isUsable())
            {
                return *particles_->levelGhostParticlesNew;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }



        NO_DISCARD field_type const& density() const
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


        NO_DISCARD field_type& density()
        {
            return const_cast<field_type&>(static_cast<const IonPopulation*>(this)->density());
        }



        NO_DISCARD VecField const& flux() const { return flux_; }
        NO_DISCARD VecField& flux() { return flux_; }

        NO_DISCARD TensorField const& momentumTensor() const { return momentumTensor_; }
        NO_DISCARD TensorField& momentumTensor() { return momentumTensor_; }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------



        struct MomentsProperty
        {
            std::string name;
            typename HybridQuantity::Scalar qty;
        };

        using MomentProperties = std::array<MomentsProperty, 1>;



        NO_DISCARD MomentProperties getFieldNamesAndQuantities() const
        {
            return {{{name_ + "_rho", HybridQuantity::Scalar::rho}}};
        }



        struct ParticleProperty
        {
            std::string name;
        };



        using ParticleProperties = std::array<ParticleProperty, 1>;

        NO_DISCARD ParticleProperties getParticleArrayNames() const { return {{{name_}}}; }




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



        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(flux_, momentumTensor_);
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
        field_type* rho_{nullptr};
        ParticlesPack<ParticleArray>* particles_{nullptr};
        initializer::PHAREDict const& particleInitializerInfo_;
    };
} // namespace core
} // namespace PHARE

#endif
