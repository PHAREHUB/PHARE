#ifndef PHARE_IONS_H
#define PHARE_IONS_H

#include <algorithm>
#include <functional>
#include <iterator>

#include "data/ions/ion_population/ion_population.h"
#include "data_provider.h"
#include "hybrid/hybrid_quantities.h"
#include "particle_initializers/particle_initializer_factory.h"


namespace PHARE
{
namespace core
{
    template<typename IonPopulation, typename GridLayout>
    class Ions
    {
    public:
        using field_type          = typename IonPopulation::field_type;
        using vecfield_type       = typename IonPopulation::vecfield_type;
        using particle_array_type = typename IonPopulation::particle_array_type;
        using ParticleInitializerFactoryT
            = ParticleInitializerFactory<particle_array_type, GridLayout>;

        static constexpr auto dimension = GridLayout::dimension;




        explicit Ions(PHARE::initializer::PHAREDict<dimension> dict)
            : name_{dict["name"].template to<std::string>()}
            , bulkVelocity_{name_ + "_bulkVel", HybridQuantity::Vector::V}
            , populations_{}
        {
            auto nbrPop = dict["nbrPopulations"].template to<std::size_t>();
            populations_.reserve(nbrPop);

            for (uint32 ipop = 0; ipop < nbrPop; ++ipop)
            {
                auto& pop = dict["pop" + std::to_string(ipop)];
                populations_.push_back(IonPopulation{name_, pop});
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

                auto& popDensity = pop.density();
                std::transform(std::begin(*rho_), std::end(*rho_), std::begin(popDensity),
                               std::begin(*rho_), std::plus<typename field_type::type>{});
            }
        }


        // TODO 3347 compute ion bulk velocity

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
