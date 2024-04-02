#ifndef PHARE_ELECTRON_POPULATION_HPP
#define PHARE_ELECTRON_POPULATION_HPP

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <array>


#include "core/def.hpp"
#include "core/pic/pic_quantities.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"

// Only works for a single population of electrons, but could be modified for many.
namespace PHARE::core
{
    template<typename ParticleArray, typename VecField, typename TensorField, typename GridLayout>
    class ElectronPopulation
    {
    public:
        using field_type                       = typename VecField::field_type;
        static constexpr std::size_t dimension = VecField::dimension;
        using particle_array_type              = ParticleArray;
        using particle_resource_type           = ParticlesPack<ParticleArray>;
        using vecfield_type                    = VecField;
        using tensorfield_type                 = TensorField;
        double Me_ov_Mp = 1/1836.152673 // mass of an electron normalized by proton mass
        std::string defName_ = PICelectrons; 

        ElectronPopulation()
            : name_{defName_}
            , mass_{Me_ov_Mp}
            , flux_{name_ + "_flux", PICQuantity::Vector::Ve} // TODO: should this is kept (for coarsening)?
        {
        } 


        NO_DISCARD double mass() const { return mass_; }

        NO_DISCARD std::string const& name() const { return name_; }


        initializer::PHAREDict initElectrons()
        {
        initializer::PHAREDict dict;

        dict["name"]    = std::string{"maxwellian"};
        dict["density"] = 0;

        dict["bulk_velocity_x"] = 0;
        dict["bulk_velocity_y"] = 0;
        dict["bulk_velocity_z"] = 0;

        dict["thermal_velocity_x"] = 0;
        dict["thermal_velocity_y"] = 0;
        dict["thermal_velocity_z"] = 0;

        dict["nbr_part_per_cell"] = 100;
        dict["charge"]            = 1.;
        dict["basis"]             = std::string{"cartesian"};
        
        return dict;
        }

        static auto init_dict_ = initElectrons();

        NO_DISCARD auto const& particleInitializerInfo() const { return init_dict_; }


        NO_DISCARD bool isUsable() const
        {
            return particles_ != nullptr && rho_ != nullptr && flux_.isUsable();
        }


        NO_DISCARD bool isSettable() const
        {
            return particles_ == nullptr && rho_ == nullptr && flux_.isSettable();
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
                static_cast<ElectronPopulation const*>(this)->domainParticles());
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
                static_cast<ElectronPopulation const*>(this)->patchGhostParticles());
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
                static_cast<ElectronPopulation const*>(this)->levelGhostParticles());
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
            return const_cast<field_type&>(static_cast<const ElectronPopulation*>(this)->density());
        }



        NO_DISCARD VecField const& flux() const { return flux_; }
        NO_DISCARD VecField& flux() { return flux_; }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------

        struct MomentsProperty
        {
            std::string name;
            typename PICQuantity::Scalar qty;
        };

        using MomentProperties = std::array<MomentsProperty, 1>;



        NO_DISCARD MomentProperties getFieldNamesAndQuantities() const
        {
            return {{{name_ + "_rho", PICQuantity::Scalar::rho}}};
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
                rho_ = field;
            else
                throw std::runtime_error("Error - invalid density buffer name");
        }



        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return flux_;
        }
        


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------



        NO_DISCARD std::string to_str()
        {
            std::stringstream ss;
            ss << "Electron Population\n";
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
    };
} // namespace PHARE::core

#endif
