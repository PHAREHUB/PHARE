#ifndef PHARE_SPLIT_H
#define PHARE_SPLIT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "data/grid/gridlayout.h"
#include "data/particles/particle.h"
#include "utilities/box/box.h"
#include "utilities/types.h"

namespace PHARE
{
#define ISNOTTABULATED(dim, RF) (dim != 1 || RF != 2)

template<std::size_t dimension, std::size_t interpOrder>
class Split
{
private:
    uint32 refinedParticlesNbr_;
    std::vector<float> weights_;
    std::vector<uint32> iCellsX_;
    std::vector<float> deltasX_;
    uint32 refinementFactor_;


    // dimension = 1, refinement factor = 2, nbrOfBabies = 2
    constexpr static std::array<float, 3> tabD1RF2N02Weight_ = {{0.5, 0.5, 0.5}};
    constexpr static std::array<float, 3> tabD1RF2N02Delta_  = {{0.277f, 0.332f, 0.376f}};


    // dimension = 1, refinement factor = 2, nbrOfBabies = 3
    constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N03Weight_
        = {{{{0.5f, 0.468f, 0.474f}}, {{0.25f, 0.266f, 0.263f}}}};
    constexpr static std::array<float, 3> tabD1RF2N03Delta_ = {{0.5f, 0.556f, 0.638f}};


    // dimension = 1, refinement factor = 2, nbrOfBabies = 4
    constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N04Weight_
        = {{{{0.0f, 0.125f, 0.135f}}, {{0.0f, 0.375f, 0.365f}}}};
    constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N04Delta_
        = {{{{0.0f, 0.75f, 0.833f}}, {{0.0f, 0.25f, 0.272f}}}};


    // dimension = 1, refinement factor = 2, nbrOfBabies = 5
    constexpr static std::array<std::array<float, 3>, 3> tabD1RF2N05Weight_
        = {{{{0.0f, 0.0f, 0.375f}}, {{0.0f, 0.0f, 0.0625f}}, {{0.0f, 0.0f, 0.25f}}}};
    constexpr static std::array<std::array<float, 3>, 2> tabD1RF2N05Delta_
        = {{{{0.0f, 0.0f, 1.f}}, {{0.0f, 0.0f, 0.5f}}}};

    const std::array<std::vector<int>, 3> tabNbrOfBabies_ = {{{2, 3}, {2, 3, 4}, {2, 3, 4, 5}}};




public:
    Split(uint32 refineFactor, uint32 refinedParticlesNbr)
        : refinedParticlesNbr_{refinedParticlesNbr}
        , refinementFactor_{refineFactor}
    {
        deltasX_.assign(refinedParticlesNbr, 0);
        weights_.assign(refinedParticlesNbr, 0);

        std::vector<int> const& acceptedNbrOfBabies = tabNbrOfBabies_[interpOrder - 1];

        if (ISNOTTABULATED(dimension, refineFactor))
        {
            std::cout << "dimension and/or refinement factor not tabulated" << std::endl;
        }

        else
        {
            if (std::find(acceptedNbrOfBabies.begin(), acceptedNbrOfBabies.end(),
                          refinedParticlesNbr)
                == acceptedNbrOfBabies.end())
            {
                std::cout << "# of babies for splitting not correct" << std::endl;
            }

            else
            {
                // weights & deltas are coming from the tabulated values
                switch (refinedParticlesNbr)
                {
                    case 2:
                        weights_[0] = tabD1RF2N02Weight_[interpOrder - 1];
                        weights_[1] = tabD1RF2N02Weight_[interpOrder - 1];

                        deltasX_[0] = -tabD1RF2N02Delta_[interpOrder - 1];
                        deltasX_[1] = +tabD1RF2N02Delta_[interpOrder - 1];
                        break;

                    case 3:
                        weights_[0] = tabD1RF2N03Weight_[0][interpOrder - 1];
                        weights_[1] = tabD1RF2N03Weight_[1][interpOrder - 1];
                        weights_[2] = tabD1RF2N03Weight_[1][interpOrder - 1];

                        deltasX_[0] = 0.0;
                        deltasX_[1] = -tabD1RF2N03Delta_[interpOrder - 1];
                        deltasX_[2] = +tabD1RF2N03Delta_[interpOrder - 1];
                        break;

                    case 4:
                        weights_[0] = tabD1RF2N04Weight_[0][interpOrder - 1];
                        weights_[1] = tabD1RF2N04Weight_[0][interpOrder - 1];
                        weights_[2] = tabD1RF2N04Weight_[1][interpOrder - 1];
                        weights_[3] = tabD1RF2N04Weight_[1][interpOrder - 1];

                        deltasX_[0] = -tabD1RF2N04Delta_[0][interpOrder - 1];
                        deltasX_[1] = +tabD1RF2N04Delta_[0][interpOrder - 1];
                        deltasX_[2] = -tabD1RF2N04Delta_[1][interpOrder - 1];
                        deltasX_[3] = +tabD1RF2N04Delta_[1][interpOrder - 1];
                        break;

                    case 5:
                        weights_[0] = tabD1RF2N05Weight_[0][interpOrder - 1];
                        weights_[1] = tabD1RF2N05Weight_[1][interpOrder - 1];
                        weights_[2] = tabD1RF2N05Weight_[1][interpOrder - 1];
                        weights_[3] = tabD1RF2N05Weight_[2][interpOrder - 1];
                        weights_[4] = tabD1RF2N05Weight_[2][interpOrder - 1];

                        deltasX_[0] = 0.0;
                        deltasX_[1] = -tabD1RF2N05Delta_[0][interpOrder - 1];
                        deltasX_[2] = +tabD1RF2N05Delta_[0][interpOrder - 1];
                        deltasX_[3] = -tabD1RF2N05Delta_[1][interpOrder - 1];
                        deltasX_[4] = +tabD1RF2N05Delta_[1][interpOrder - 1];
                        break;

                    default: std::cout << "ta mere en short !" << std::endl;
                }
            }
        }
    }

    ~Split() = default;


    template<typename GridLayoutT>
    inline void operator()(Particle<dimension> const& coarsePartOnRefinedGrid,
                           std::vector<Particle<dimension>>& refinedParticles) const
    {
        for (uint32 refinedParticleIndex = 0; refinedParticleIndex < refinedParticlesNbr_;
             ++refinedParticleIndex)
        {
            if constexpr (dimension == 1)
            {
                // the values for icell & delta are only working for 1 dim...
                float weight = coarsePartOnRefinedGrid.weight * weights_[refinedParticleIndex];
                int32 icell  = coarsePartOnRefinedGrid.icell[0];
                float delta  = coarsePartOnRefinedGrid.delta[0]
                              + deltasX_[refinedParticleIndex] * refinementFactor_;

                // weights & deltas are the only known values for the babies.
                // so the icell values of each baby needs to be calculated
                float integra = std::floor(delta);
                delta -= integra;
                icell += static_cast<int32>(integra);

                refinedParticles.emplace_back(weight, coarsePartOnRefinedGrid.charge, {{icell}},
                                              {{delta}}, coarsePartOnRefinedGrid.v);
            }
            else if constexpr (dimension != 1)
            {
                static_assert("Only 1D is supported for split at the moment");
            }
        }
    }
};




} // namespace PHARE
#endif // endif SPLIT_H
