#ifndef PHARE_IONS_HPP
#define PHARE_IONS_HPP

#include <algorithm>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>
#include <cmath>


#include "core/def.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "initializer/data_provider.hpp"
#include "particle_initializers/particle_initializer_factory.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

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
        using Float               = typename field_type::type;
        using particle_array_type = typename IonPopulation::particle_array_type;
        using ParticleInitializerFactoryT
            = ParticleInitializerFactory<particle_array_type, GridLayout>;
        using gridlayout_type           = GridLayout;
        static constexpr auto dimension = GridLayout::dimension;




        explicit Ions(PHARE::initializer::PHAREDict const& dict)
            : bulkVelocity_{"bulkVel", HybridQuantity::Vector::V}
            , populations_{generate(
                  [&dict](auto ipop) { //
                      return IonPopulation{dict["pop" + std::to_string(ipop)]};
                  },
                  dict["nbrPopulations"].template to<std::size_t>())}
            , sameMasses_{allSameMass_()}
        {
        }


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
        NO_DISCARD std::string massDensityName() const { return "massDensity"; }


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
        void computeMassDensity()
        {
            massDensity_->zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto& popDensity = pop.density();
                std::transform(
                    std::begin(*massDensity_), std::end(*massDensity_), std::begin(popDensity),
                    std::begin(*massDensity_),
                    [&pop](auto const& n, auto const& pop_n) { return n + pop_n * pop.mass(); });
            }
        }


        void computeBulkVelocity()
        {
            // the bulk velocity is sum(pop_mass * pop_flux) / sum(pop_mass * pop_density)
            // if all populations have the same mass, this is equivalent to sum(pop_flux) /
            // sum(pop_density) sum(pop_density) is rho_ and already known by the time we get here.
            // sum(pop_mass * pop_flux) is massDensity_ and is computed by computeMassDensity() if
            // needed
            if (!sameMasses_)
                computeMassDensity();
            auto const& density = (sameMasses_) ? rho_ : massDensity_;

            bulkVelocity_.zero();
            auto& vx = bulkVelocity_.getComponent(Component::X);
            auto& vy = bulkVelocity_.getComponent(Component::Y);
            auto& vz = bulkVelocity_.getComponent(Component::Z);

            for (auto& pop : populations_)
            {
                // account for mass only if populations have different masses
                std::function<Float(Float, Float)> plus = std::plus<Float>{};
                std::function<Float(Float, Float)> plusMass
                    = [&pop](Float const& v, Float const& f) { return v + f * pop.mass(); };

                auto const& flux    = pop.flux();
                auto&& [fx, fy, fz] = flux();


                std::transform(std::begin(vx), std::end(vx), std::begin(fx), std::begin(vx),
                               sameMasses_ ? plus : plusMass);
                std::transform(std::begin(vy), std::end(vy), std::begin(fy), std::begin(vy),
                               sameMasses_ ? plus : plusMass);
                std::transform(std::begin(vz), std::end(vz), std::begin(fz), std::begin(vz),
                               sameMasses_ ? plus : plusMass);
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


        // in the following isUsable and isSettable the massDensity_ is not checked
        // because it is for internal use only so no object will ever need to access it.
        NO_DISCARD bool isUsable() const
        {
            bool usable = rho_ != nullptr && bulkVelocity_.isUsable();
            for (auto const& pop : populations_)
            {
                usable = usable && pop.isUsable();
            }
            return usable;
        }



        NO_DISCARD bool isSettable() const
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

        NO_DISCARD MomentProperties getFieldNamesAndQuantities() const
        {
            return {{{densityName(), HybridQuantity::Scalar::rho},
                     {massDensityName(), HybridQuantity::Scalar::rho}}};
        }



        void setBuffer(std::string const& bufferName, field_type* field)
        {
            if (bufferName == densityName())
            {
                rho_ = field;
            }
            else if (bufferName == massDensityName())
            {
                massDensity_ = field;
            }
            else
            {
                throw std::runtime_error("Error - invalid density buffer name : " + bufferName);
            }
        }



        NO_DISCARD std::vector<IonPopulation>& getRunTimeResourcesUserList()
        {
            return populations_;
        }

        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(bulkVelocity_);
        }



        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------


        NO_DISCARD std::string to_str()
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
        bool allSameMass_() const
        {
            return all(populations_, [this](auto const& pop) { // arbitrary small diff
                return float_equals(pop.mass(), populations_.front().mass(), /*abs_tol=*/1e-10);
            });
        }




        field_type* rho_{nullptr};
        field_type* massDensity_{nullptr};
        vecfield_type bulkVelocity_;
        std::vector<IonPopulation> populations_;

        // note this is set at construction when all populations are added
        // in the future if some populations are dynamically created during the simulation
        // this should be updated accordingly
        bool sameMasses_{false};
    };
} // namespace core
} // namespace PHARE

#endif
