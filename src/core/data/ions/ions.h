#ifndef PHARE_IONS_H
#define PHARE_IONS_H

#include <algorithm>
#include <functional>
#include <iterator>

#include "data/ions/ion_population/ion_population.h"
#include "hybrid/hybrid_quantities.h"
#include "ion_initializer.h"

namespace PHARE
{
namespace core
{
    template<typename IonPopulation, typename GridLayout>
    class Ions
    {
    public:
        using field_type    = typename IonPopulation::field_type;
        using vecfield_type = typename IonPopulation::vecfield_type;
        using ions_initializer_type
            = IonsInitializer<typename IonPopulation::particle_array_type, GridLayout>;

        explicit Ions(ions_initializer_type initializer)
            : name_{std::move(initializer.name)}
            , bulkVelocity_{name_ + "_bulkVel", HybridQuantity::Vector::V}
            , populations_{}
        {
            // TODO IonPopulation constructor will need to take a ParticleInitializer
            // from the vector in the initializer
            populations_.reserve(initializer.nbrPopulations);
            for (uint32 ipop = 0; ipop < initializer.nbrPopulations; ++ipop)
            {
                populations_.push_back(
                    IonPopulation{name_ + "_" + initializer.names[ipop], initializer.masses[ipop],
                                  std::move(initializer.particleInitializers[ipop])});
            }
        }


        /*
         *         explicit Ions(PHARE::initializer::PHAREDict<1>& dict)
            : name_{dict["name"].to<std::string>()}
            , bulkVelocity_{name_ + "_bulkVel", HybridQuantity::Vector::V}
            , populations_{}
        {
            // TODO IonPopulation constructor will need to take a ParticleInitializer
            // from the vector in the initializer
            populations_.reserve(dict["nbrPopulation"].to<std::size_t>());
            // initializer.nbrPopulations);
            for (uint32 ipop = 0; ipop < populations_.size(); ++ipop)
            {
                auto& pop    = dict["pop" + std::to_string(ipop)];
                auto popName = pop["name"].to<std::string>();
                auto mass    = pop["mass"].to<double>();
                // auto paticleInitializer = std::move();

                // TODO use ParticleInitializerFactory with subdict info
                populations_.push_back(name_ + "_" + popName, mass, nullptr);
            }
        }
         */

        field_type const& density() const
        {
            if (isUsable())
            {
                return *rho_;
            }
            else
            {
                throw std::runtime_error("Error - cannot access density data");
            }
        }



        field_type& density()
        {
            if (isUsable())
            {
                return *rho_;
            }
            else
            {
                throw std::runtime_error("Error - cannot access density data");
            }
        }




        void loadParticles(GridLayout const& layout)
        {
            for (auto& pop : populations_)
            {
                pop.loadParticles(layout);
            }
        }


        vecfield_type const& velocity() const { return bulkVelocity_; }

        vecfield_type& velocity() { return bulkVelocity_; }

        std::string densityName() const { return name_ + "_rho"; }


        void computeDensity()
        {
            rho_->zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto const& popDensity = pop.density();
                std::transform(std::begin(rho_), std::end(rho_), std::begin(rho_),
                               std::plus<typename field_type::type>{});
            }
        }



        auto begin() { return std::begin(populations_); }
        auto end() { return std::end(populations_); }

        auto begin() const { return std::begin(populations_); }
        auto end() const { return std::end(populations_); }


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
            return {{{densityName(), HybridQuantity::Scalar::rho}}};
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
} // namespace core
} // namespace PHARE

#endif
