#ifndef PHARE_TYPES_H
#define PHARE_TYPES_H

#include "include.h"

namespace PHARE
{
template<std::size_t dimension_, std::size_t interp_order_, std::size_t nbRefinedPart_>
struct PHARE_Types
{
    static auto constexpr dimension     = dimension_;
    static auto constexpr interp_order  = interp_order_;
    static auto constexpr nbRefinedPart = nbRefinedPart_;

    using Array_t
        = decltype(PHARE::core::makeNdArray<dimension>(std::array<std::uint32_t, dimension>{}));
    using VecField_t      = PHARE::core::VecField<Array_t, PHARE::core::HybridQuantity>;
    using Field_t         = PHARE::core::Field<Array_t, PHARE::core::HybridQuantity::Scalar>;
    using Electromag_t    = PHARE::core::Electromag<VecField_t>;
    using YeeLayout_t     = PHARE::core::GridLayoutImplYee<dimension, interp_order>;
    using GridLayout_t    = PHARE::core::GridLayout<YeeLayout_t>;
    using ParticleArray_t = PHARE::core::ParticleArray<dimension>;
    using MaxwellianParticleInitializer_t
        = PHARE::core::MaxwellianParticleInitializer<ParticleArray_t, GridLayout_t>;
    using IonPopulation_t = PHARE::core::IonPopulation<ParticleArray_t, VecField_t, GridLayout_t>;
    using Ions_t          = PHARE::core::Ions<IonPopulation_t, GridLayout_t>;
    using Electrons_t     = PHARE::core::Electrons<Ions_t>;
    using IPhysicalModel  = PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>;
    using HybridModel_t   = PHARE::solver::HybridModel<GridLayout_t, Electromag_t, Ions_t,
                                                     Electrons_t, PHARE::amr::SAMRAI_Types>;
    using MHDModel_t  = PHARE::solver::MHDModel<GridLayout_t, VecField_t, PHARE::amr::SAMRAI_Types>;
    using SolverPPC_t = PHARE::solver::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;
    using SolverMHD_t = PHARE::solver::SolverMHD<MHDModel_t, PHARE::amr::SAMRAI_Types>;
    using LevelInitializerFactory_t = PHARE::solver::LevelInitializerFactory<HybridModel_t>;
    using hierarchy_t               = typename PHARE::amr::SAMRAI_Types::hierarchy_t;
    using ModelView                 = PHARE::diagnostic::ModelView<hierarchy_t, HybridModel_t>;
    using DiagnosticWriter          = PHARE::diagnostic::h5::Writer<ModelView>;

    using ParticleInitializerFactory
        = PHARE::core::ParticleInitializerFactory<ParticleArray_t, GridLayout_t>;
    using Splitter         = PHARE::amr::Split<dimension, interp_order, nbRefinedPart>;
    using RefinementParams = PHARE::amr::RefinementParams<Splitter>;
    using MessengerFactory
        = PHARE::amr::MessengerFactory<MHDModel_t, HybridModel_t, IPhysicalModel, RefinementParams>;

    using MultiPhysicsIntegrator
        = PHARE::solver::MultiPhysicsIntegrator<MessengerFactory, LevelInitializerFactory_t,
                                                PHARE::amr::SAMRAI_Types>;
};

} // namespace PHARE

#endif // PHARE_TYPES_H
