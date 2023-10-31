#ifndef PHARE_ION_POPULATION_HPP
#define PHARE_ION_POPULATION_HPP

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>


#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/ions/particle_initializers/particle_initializer.hpp"
#include "initializer/data_provider.hpp"
#include "particle_pack.hpp"
#include "core/utilities/algorithm.hpp"

namespace PHARE
{
namespace core
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


        IonPopulation(initializer::PHAREDict const& initializer)
            : name_{initializer["name"].template to<std::string>()}
            , mass_{initializer["mass"].template to<double>()}
            , flux_{name_ + "_flux", HybridQuantity::Vector::V}
            , particleInitializerInfo_{initializer["particle_initializer"]}
        {
        }


        [[nodiscard]] double mass() const { return mass_; }

        [[nodiscard]] std::string const& name() const { return name_; }


        [[nodiscard]] auto const& particleInitializerInfo() const
        {
            return particleInitializerInfo_;
        }



        [[nodiscard]] bool isUsable() const
        {
            return particles_ != nullptr && rho_ != nullptr && flux_.isUsable();
        }


        [[nodiscard]] bool isSettable() const
        {
            return particles_ == nullptr && rho_ == nullptr && flux_.isSettable();
        }




        [[nodiscard]] auto nbrParticles() const
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




        [[nodiscard]] auto& domainParticles() const
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

        [[nodiscard]] auto& domainParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->domainParticles());
        }



        [[nodiscard]] auto& patchGhostParticles() const
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

        [[nodiscard]] auto& patchGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->patchGhostParticles());
        }


        [[nodiscard]] auto& levelGhostParticles() const
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

        [[nodiscard]] auto& levelGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->levelGhostParticles());
        }




        [[nodiscard]] ParticleArray& levelGhostParticlesOld()
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



        [[nodiscard]] ParticleArray& levelGhostParticlesNew()
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



        [[nodiscard]] field_type const& density() const
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


        [[nodiscard]] field_type& density()
        {
            return const_cast<field_type&>(static_cast<const IonPopulation*>(this)->density());
        }



        [[nodiscard]] VecField const& flux() const { return flux_; }


        [[nodiscard]] VecField& flux() { return flux_; }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------



        struct MomentsProperty
        {
            std::string name;
            typename HybridQuantity::Scalar qty;
        };

        using MomentProperties = std::vector<MomentsProperty>;



        [[nodiscard]] MomentProperties getFieldNamesAndQuantities() const
        {
            return {{{name_ + "_rho", HybridQuantity::Scalar::rho}}};
        }



        struct ParticleProperty
        {
            std::string name;
        };



        using ParticleProperties = std::vector<ParticleProperty>;

        [[nodiscard]] ParticleProperties getParticleArrayNames() const { return {{{name_}}}; }




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



        [[nodiscard]] auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(flux_);
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------



        [[nodiscard]] std::string to_str()
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
        field_type* rho_{nullptr};
        ParticlesPack<ParticleArray>* particles_{nullptr};
        initializer::PHAREDict const& particleInitializerInfo_;
    };
} // namespace core
} // namespace PHARE

#endif
