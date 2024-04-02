#ifndef PHARE_PIC_ELECTRONS_HPP
#define PHARE_PIC_ELECTRONS_HPP

#include <algorithm>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <array>


#include "core/def.hpp"
#include "core/pic/pic_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "initializer/data_provider.hpp"
#include "particle_initializers/particle_initializer_factory.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/data/fermions/electron_population.hpp"

namespace PHARE::core
{
    template<typename ElectronPopulation, typename GridLayout>
    class PICElectrons
    {
    public:
        using field_type          = typename ElectronPopulation::field_type;
        using vecfield_type       = typename ElectronPopulation::vecfield_type;
        using Float               = typename field_type::type;
        using tensorfield_type    = typename ElectronPopulation::tensorfield_type;
        using particle_array_type = typename ElectronPopulation::particle_array_type;
        using ParticleInitializerFactoryT
            = ParticleInitializerFactory<particle_array_type, GridLayout>;
        using gridlayout_type           = GridLayout;
        static constexpr auto dimension = GridLayout::dimension;



        explicit PICElectrons()
            : bulkVelocity_{"ElectronBulkVel", PICQuantity::Vector::Ve}
            , populations_{} // TODO check how to initialize without using PHAREdict
        {
        }

        auto const& populations_ = ElectronPopulation{}; // Placeholder till I find smth better


        NO_DISCARD auto nbrPopulations() const { return populations_.size(); }


        NO_DISCARD field_type const& density() const
        {
            if (isUsable())
                return *rho_;
            else
                throw std::runtime_error("Error - cannot access density data");
        }



        NO_DISCARD field_type& density()
        {
            if (isUsable())
                return *rho_;
            else
                throw std::runtime_error("Error - cannot access density data");
        }



        NO_DISCARD vecfield_type const& velocity() const { return bulkVelocity_; }

        NO_DISCARD vecfield_type& velocity() { return bulkVelocity_; }

        NO_DISCARD std::string densityName() const { return "rho"; }

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
                               std::begin(*rho_), std::plus<Float>{});
            }
        }



        void computeBulkVelocity()
        {
            // the bulk velocity is sum(pop_mass * pop_flux) / sum(pop_mass * pop_density)
            // Since all populations have the same mass, this is equivalent to sum(pop_flux) /
            // sum(pop_density) sum(pop_density) which is rho_ and already known by the time we get here.
 
            auto const& density = rho_ ;

            bulkVelocity_.zero();
            auto& vx = bulkVelocity_.getComponent(Component::X);
            auto& vy = bulkVelocity_.getComponent(Component::Y);
            auto& vz = bulkVelocity_.getComponent(Component::Z);

            for (auto& pop : populations_)
            {
                std::function<Float(Float, Float)> plus = std::plus<Float>{};

                auto const& flux    = pop.flux();
                auto&& [fx, fy, fz] = flux();


                std::transform(std::begin(vx), std::end(vx), std::begin(fx), std::begin(vx),
                               plus);
                std::transform(std::begin(vy), std::end(vy), std::begin(fy), std::begin(vy),
                               plus);
                std::transform(std::begin(vz), std::end(vz), std::begin(fz), std::begin(vz),
                               plus);
            }


            std::transform(std::begin(vx), std::end(vx), std::begin(*density), std::begin(vx),
                           std::divides<Float>{});
            std::transform(std::begin(vy), std::end(vy), std::begin(*density), std::begin(vy),
                           std::divides<Float>{});
            std::transform(std::begin(vz), std::end(vz), std::begin(*density), std::begin(vz),
                           std::divides<Float>{});
        }



        NO_DISCARD auto begin() { return std::begin(populations_); }
        NO_DISCARD auto end() { return std::end(populations_); }
        NO_DISCARD auto begin() const { return std::begin(populations_); }
        NO_DISCARD auto end() const { return std::end(populations_); }


        NO_DISCARD bool isUsable() const
        {
            bool usable
                = rho_ != nullptr && bulkVelocity_.isUsable();

            for (auto const& pop : populations_)
            {
                usable = usable && pop.isUsable();
            }
            return usable;
        }



        NO_DISCARD bool isSettable() const
        {
            bool settable
                = rho_ == nullptr && bulkVelocity_.isSettable()

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
            typename PICQuantity::Scalar qty;
        };

        using MomentProperties = std::vector<MomentsProperty>;

        NO_DISCARD MomentProperties getFieldNamesAndQuantities() const
        {
            return {{{densityName(), PICQuantity::Scalar::rho}}};
        }



        void setBuffer(std::string const& bufferName, field_type* field)
        {
            if (bufferName == densityName())
                rho_ = field;
            else
                throw std::runtime_error("Error - invalid density buffer name : " + bufferName);
        }



        NO_DISCARD std::vector<ElectronPopulation>& getRunTimeResourcesUserList()
        {
            return populations_;
        }

        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return bulkVelocity_;
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------


        NO_DISCARD std::string to_str()
        {
            std::stringstream ss;
            ss << "PICElectrons\n";
            ss << "------------------------------------\n";
            ss << "number of populations  : " << nbrPopulations() << "\n";
            for (auto& pop : populations_)
                ss << core::to_str(pop);
            return ss.str();
        }



        field_type* rho_{nullptr};
        vecfield_type bulkVelocity_;
        std::vector<ElectronPopulation> populations_;
    };
} // namespace PHARE::core

#endif
