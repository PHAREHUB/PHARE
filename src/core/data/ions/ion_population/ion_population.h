#ifndef PHARE_ION_POPULATION_H
#define PHARE_ION_POPULATION_H

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>


#include "core/hybrid/hybrid_quantities.h"
#include "core/data/ions/particle_initializers/particle_initializer.h"
#include "initializer/data_provider.h"
#include "particle_pack.h"
#include "core/utilities/algorithm.h"
#include "core/utilities/types.h"

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


        double mass() const { return mass_; }

        std::string const& name() const { return name_; }


        auto const& particleInitializerInfo() const { return particleInitializerInfo_; }



        bool isUsable() const
        {
            return particles_ != nullptr && rho_ != nullptr && flux_.isUsable();
        }


        bool isSettable() const
        {
            return particles_ == nullptr && rho_ == nullptr && flux_.isSettable();
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




        auto& domainParticles() const
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

        auto& domainParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->domainParticles());
        }



        auto& patchGhostParticles() const
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

        auto& patchGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->patchGhostParticles());
        }


        auto& levelGhostParticles() const
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

        auto& levelGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<IonPopulation const*>(this)->levelGhostParticles());
        }



        ParticleArray& levelGhostParticlesOld()
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



        ParticleArray& levelGhostParticlesNew()
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



        VecField const& flux() const { return flux_; }


        VecField& flux() { return flux_; }



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



        std::string to_str()
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

namespace PHARE::core
{
template<typename ParticleArray, typename VecField, typename GridLayout_>
struct IonPopulationView
{
    using This = IonPopulationView<ParticleArray, VecField, GridLayout_>;

    using particle_array_type = ParticleArray;
    using VecFieldView_t      = typename VecField::view_t;
    using Field_t             = typename VecField::field_type;
    using FieldView_t         = typename Field_t::view_t;

    FieldView_t density;
    VecFieldView_t flux;
    ParticleArray* domain;
    ParticleArray* patch_ghost;
    ParticleArray* level_ghost;
    double const mass;


    auto& domainParticles() { return *domain; }


    std::shared_ptr<This> static make_shared(
        IonPopulation<ParticleArray, VecField, GridLayout_>& pop)
    {
        return std::make_shared<aggregate_adapter<This>>(
            pop.density().view(), pop.flux().view(), &pop.domainParticles(),
            &pop.patchGhostParticles(), &pop.levelGhostParticles(), pop.mass());
    }

    template<typename Ions>
    auto static make_shared(Ions& ions)
    {
        return generate([&](auto& pop) { return make_shared(pop); }, ions);
    }
};
} // namespace PHARE::core

#endif
