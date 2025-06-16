#ifndef PHARE_IONS_HPP
#define PHARE_IONS_HPP

#include <algorithm>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <array>


#include "core/def.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "initializer/data_provider.hpp"
#include "particle_initializers/particle_initializer_factory.hpp"
#include "core/utilities/algorithm.hpp"

namespace PHARE
{
namespace core
{
    template<typename IonPopulation, typename GridLayout>
    class Ions
    {
    public:
        using value_type                = IonPopulation;
        using field_type                = typename IonPopulation::field_type;
        using vecfield_type             = typename IonPopulation::vecfield_type;
        using Float                     = typename field_type::type;
        using tensorfield_type          = typename IonPopulation::tensorfield_type;
        using particle_array_type       = typename IonPopulation::particle_array_type;
        using gridlayout_type           = GridLayout;
        static constexpr auto dimension = GridLayout::dimension;


        Ions(Ions const&) = default;
        Ions(Ions&&)      = default;


        explicit Ions(PHARE::initializer::PHAREDict const& dict)
            : massDensity_{massDensityName(), HybridQuantity::Scalar::rho}
            , chargeDensity_{chargeDensityName(), HybridQuantity::Scalar::rho}
            , bulkVelocity_{"bulkVel", HybridQuantity::Vector::V}
            , populations_{generate(
                  [&dict](auto ipop) { //
                      return IonPopulation{dict["pop" + std::to_string(ipop)]};
                  },
                  dict["nbrPopulations"].template to<std::size_t>())}
            , momentumTensor_{"momentumTensor", HybridQuantity::Tensor::M}
        {
        }


        NO_DISCARD auto nbrPopulations() const { return populations_.size(); }
        NO_DISCARD auto size() const { return nbrPopulations(); }

        NO_DISCARD field_type const& massDensity() const { return massDensity_; }
        NO_DISCARD field_type const& massDensity() { return massDensity_; }

        NO_DISCARD field_type const& chargeDensity() const { return chargeDensity_; }
        NO_DISCARD field_type& chargeDensity() { return chargeDensity_; }

        NO_DISCARD vecfield_type const& velocity() const { return bulkVelocity_; }
        NO_DISCARD vecfield_type& velocity() { return bulkVelocity_; }

        NO_DISCARD std::string static chargeDensityName() { return "chargeDensity"; }
        NO_DISCARD std::string static massDensityName() { return "massDensity"; }

        tensorfield_type const& momentumTensor() const { return momentumTensor_; }
        tensorfield_type& momentumTensor() { return momentumTensor_; }

        void computeChargeDensity()
        {
            chargeDensity_.zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto& popDensity = pop.chargeDensity();
                std::transform(std::begin(chargeDensity_), std::end(chargeDensity_),
                               std::begin(popDensity), std::begin(chargeDensity_),
                               std::plus<Float>{});
            }
        }

        void computeMassDensity()
        {
            massDensity_.zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto& popDensity = pop.particleDensity();
                std::transform(
                    std::begin(massDensity_), std::end(massDensity_), std::begin(popDensity),
                    std::begin(massDensity_),
                    [&pop](auto const& n, auto const& pop_n) { return n + pop_n * pop.mass(); });
            }
        }


        void computeBulkVelocity()
        {
            computeMassDensity();

            bulkVelocity_.zero();
            auto& vx = bulkVelocity_.getComponent(Component::X);
            auto& vy = bulkVelocity_.getComponent(Component::Y);
            auto& vz = bulkVelocity_.getComponent(Component::Z);

            for (auto& pop : populations_)
            {
                // account for mass only if populations have different masses
                std::function<Float(Float, Float)> plusMass
                    = [&pop](Float const& v, Float const& f) { return v + f * pop.mass(); };

                auto const& flux    = pop.flux();
                auto&& [fx, fy, fz] = flux();

                std::transform(std::begin(vx), std::end(vx), std::begin(fx), std::begin(vx),
                               plusMass);
                std::transform(std::begin(vy), std::end(vy), std::begin(fy), std::begin(vy),
                               plusMass);
                std::transform(std::begin(vz), std::end(vz), std::begin(fz), std::begin(vz),
                               plusMass);
            }


            std::transform(std::begin(vx), std::end(vx), std::begin(massDensity_), std::begin(vx),
                           std::divides<Float>{});
            std::transform(std::begin(vy), std::end(vy), std::begin(massDensity_), std::begin(vy),
                           std::divides<Float>{});
            std::transform(std::begin(vz), std::end(vz), std::begin(massDensity_), std::begin(vz),
                           std::divides<Float>{});
        }


        void computeFullMomentumTensor()
        {
            momentumTensor_.zero();
            auto& mom = momentumTensor_;

            for (auto& pop : populations_)
            {
                auto& p_mom = pop.momentumTensor();
                for (auto p_mij = p_mom.begin(), mij = mom.begin(); p_mij != p_mom.end();
                     ++p_mij, ++mij)
                {
                    std::transform(std::begin(*mij), std::end(*mij), std::begin(*p_mij),
                                   std::begin(*mij), std::plus<typename field_type::type>{});
                }
            }
        }


        NO_DISCARD auto begin() { return std::begin(populations_); }
        NO_DISCARD auto end() { return std::end(populations_); }
        NO_DISCARD auto begin() const { return std::begin(populations_); }
        NO_DISCARD auto end() const { return std::end(populations_); }

        NO_DISCARD auto& population(std::size_t const i)
        {
            if (i >= populations_.size())
                throw std::out_of_range("Ions population index out of range");
            return populations_[i];
        }

        NO_DISCARD auto const& population(std::size_t const i) const
        {
            if (i >= populations_.size())
                throw std::out_of_range("Ions population index out of range");
            return populations_[i];
        }

        // in the following isUsable and isSettable the massDensity_ is not checked
        // because it is for internal use only so no object will ever need to access it.
        NO_DISCARD bool isUsable() const
        {
            bool usable = chargeDensity_.isUsable() and bulkVelocity_.isUsable()
                          and momentumTensor_.isUsable() and massDensity_.isUsable();

            for (auto const& pop : populations_)
            {
                usable = usable and pop.isUsable();
            }
            return usable;
        }



        NO_DISCARD bool isSettable() const
        {
            bool settable = massDensity_.isSettable() and chargeDensity_.isSettable()
                            and bulkVelocity_.isSettable() and momentumTensor_.isSettable();

            for (auto const& pop : populations_)
            {
                settable = settable and pop.isSettable();
            }
            return settable;
        }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------



        NO_DISCARD std::vector<IonPopulation>& getRunTimeResourcesViewList()
        {
            return populations_;
        }

        NO_DISCARD auto getCompileTimeResourcesViewList()
        {
            return std::forward_as_tuple(bulkVelocity_, momentumTensor_, chargeDensity_,
                                         massDensity_);
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


        auto& operator[](std::size_t const i) const { return populations_[i]; }
        auto& operator[](std::size_t const i) { return populations_[i]; }

    private:
        field_type massDensity_;
        field_type chargeDensity_;
        vecfield_type bulkVelocity_;
        std::vector<IonPopulation> populations_;
        tensorfield_type momentumTensor_;
    };

} // namespace core
} // namespace PHARE



#endif
