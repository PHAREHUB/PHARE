#ifndef PHARE_IONS_H
#define PHARE_IONS_H

#include <algorithm>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>
#include <cmath>


#include "core/hybrid/hybrid_quantities.h"
#include "core/data/ions/ion_population/ion_population.h"
#include "core/data/vecfield/vecfield_component.h"
#include "initializer/data_provider.h"
#include "particle_initializers/particle_initializer_factory.h"
#include "core/utilities/algorithm.h"
#include "core/data/ndarray/ndarray_vector.h"

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
        using gridlayout_type           = GridLayout;
        static constexpr auto dimension = GridLayout::dimension;




        explicit Ions(PHARE::initializer::PHAREDict const& dict)
            : bulkVelocity_{"bulkVel", HybridQuantity::Vector::V}
            , populations_{}
        {
            auto nbrPop = dict["nbrPopulations"].template to<int>();

            for (int ipop = 0; ipop < nbrPop; ++ipop)
            {
                auto const& pop = dict["pop" + std::to_string(ipop)];
                populations_.emplace_back(pop);
            }
        }


        auto nbrPopulations() const { return populations_.size(); }


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

        std::string densityName() const { return "rho"; }


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

                for (std::size_t i = 5; i < rho_->size() - 5; ++i)
                {
                    assert(rho_->data()[i] > 0);
                }
            }
        }


        void computeBulkVelocity()
        {
            bulkVelocity_.zero();
            auto& vx = bulkVelocity_.getComponent(Component::X);
            auto& vy = bulkVelocity_.getComponent(Component::Y);
            auto& vz = bulkVelocity_.getComponent(Component::Z);

            for (auto& pop : populations_)
            {
                auto& flux = pop.flux();
                auto& fx   = flux.getComponent(Component::X);
                auto& fy   = flux.getComponent(Component::Y);
                auto& fz   = flux.getComponent(Component::Z);

                std::transform(std::begin(vx), std::end(vx), std::begin(fx), std::begin(vx),
                               std::plus<typename field_type::type>{});
                std::transform(std::begin(vy), std::end(vy), std::begin(fy), std::begin(vy),
                               std::plus<typename field_type::type>{});
                std::transform(std::begin(vz), std::end(vz), std::begin(fz), std::begin(vz),
                               std::plus<typename field_type::type>{});
            }

            std::transform(std::begin(vx), std::end(vx), std::begin(*rho_), std::begin(vx),
                           std::divides<typename field_type::type>{});
            std::transform(std::begin(vy), std::end(vy), std::begin(*rho_), std::begin(vy),
                           std::divides<typename field_type::type>{});
            std::transform(std::begin(vz), std::end(vz), std::begin(*rho_), std::begin(vz),
                           std::divides<typename field_type::type>{});
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
            if (bufferName == "rho")
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


        std::string to_str()
        {
            std::stringstream ss;
            ss << "Ions\n";
            ss << "------------------------------------\n";
            ss << "number of populations  : " << nbrPopulations() << "\n";
            for (auto& pop : populations_)
                ss << core::to_str(pop);
            return ss.str();
        }


    private:
        field_type* rho_{nullptr};
        vecfield_type bulkVelocity_;
        std::vector<IonPopulation> populations_; // TODO we have to name this so they are unique
                                                 // although only 1 Ions should exist.
    };



    /**
     * This routine sets NaNs on all but closest to domain ghost moment nodes
     * All nodes set to NaN are unsued and NaNs helps finding bugs if they are
     * used by mistake. The ghost node that is closest to domain is used
     * in refinement ratio 2 coarsening and should thus not be NaN. It is set
     * to a copy of the first domain node.
     *
     * Directions for start index are superfluous, as all directions are primal for these fields
     */
    template<typename Ions, typename GridLayout>
    void fixMomentGhosts(Ions& ions, GridLayout const& layout)
    {
        using Mask = NdArrayMask;

        for (auto& pop : ions)
        {
            for (auto& [id, type] : Components::componentMap)
            {
                auto& flux_comp = pop.flux().getComponent(type);

                auto phyStartIdx
                    = layout.physicalStartIndex(flux_comp.physicalQuantity(), Direction::X);

                flux_comp[Mask{0, phyStartIdx - 2}] = NAN;
                flux_comp[Mask{phyStartIdx}] >> flux_comp[Mask{phyStartIdx - 1}];
            }

            auto& density    = pop.density();
            auto phyStartIdx = layout.physicalStartIndex(density.physicalQuantity(), Direction::X);

            density[Mask{0, phyStartIdx - 2}] = NAN;
            density[Mask{phyStartIdx}] >> density[Mask{phyStartIdx - 1}];
        }
    }

} // namespace core
} // namespace PHARE

#endif
