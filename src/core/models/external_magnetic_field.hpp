#ifndef PHARE_CORE_MODELS_EXTERNAL_MAGNETIC_FIELD_HPP
#define PHARE_CORE_MODELS_EXTERNAL_MAGNETIC_FIELD_HPP

#include "core/def.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"
#include "core/mhd/mhd_quantities.hpp"

#include "initializer/data_provider.hpp"

namespace PHARE::core
{
/**
 * @brief The static background magnetic field B0 of the MHD B = B0 + B1 split.
 *
 * B0 is a known analytic function (user-supplied). It is held as a single instance by the model
 * and shared by every Runge-Kutta stage (B0 does not change between stages, so it is never
 * copied around). It carries, besides the native face-centered vector, the analytic samples of
 * B0 at every location the numerics consume it:
 *   - the constrained-transport EMF edges (atEz/atEx/atEy),
 *   - the Riemann-solver transverse face positions (atFaceX/atFaceY/atFaceZ).
 * Sampling B0 analytically at each location keeps it interpolation-free: only the perturbation
 * B1 is ever reconstructed/interpolated.
 */
template<typename VecFieldT>
class ExternalMagneticField
{
public:
    using vecfield_type = VecFieldT;
    using field_type    = typename VecFieldT::field_type;

    static constexpr auto dimension = VecFieldT::dimension;

    explicit ExternalMagneticField(PHARE::initializer::PHAREDict const& dict)
        : B0{dict["name"].template to<std::string>() + "_" + "B0", MHDQuantity::Vector::B0}
        , B0x_Ez{dict["name"].template to<std::string>() + "_" + "B0x_Ez", MHDQuantity::Scalar::Ez}
        , B0y_Ez{dict["name"].template to<std::string>() + "_" + "B0y_Ez", MHDQuantity::Scalar::Ez}
        , B0y_Ex{dict["name"].template to<std::string>() + "_" + "B0y_Ex", MHDQuantity::Scalar::Ex}
        , B0z_Ex{dict["name"].template to<std::string>() + "_" + "B0z_Ex", MHDQuantity::Scalar::Ex}
        , B0x_Ey{dict["name"].template to<std::string>() + "_" + "B0x_Ey", MHDQuantity::Scalar::Ey}
        , B0z_Ey{dict["name"].template to<std::string>() + "_" + "B0z_Ey", MHDQuantity::Scalar::Ey}
        , B0y_FaceX{dict["name"].template to<std::string>() + "_" + "B0y_FaceX",
                    MHDQuantity::Scalar::Bx}
        , B0z_FaceX{dict["name"].template to<std::string>() + "_" + "B0z_FaceX",
                    MHDQuantity::Scalar::Bx}
        , B0x_FaceY{dict["name"].template to<std::string>() + "_" + "B0x_FaceY",
                    MHDQuantity::Scalar::By}
        , B0z_FaceY{dict["name"].template to<std::string>() + "_" + "B0z_FaceY",
                    MHDQuantity::Scalar::By}
        , B0x_FaceZ{dict["name"].template to<std::string>() + "_" + "B0x_FaceZ",
                    MHDQuantity::Scalar::Bz}
        , B0y_FaceZ{dict["name"].template to<std::string>() + "_" + "B0y_FaceZ",
                    MHDQuantity::Scalar::Bz}
        , B0init_{dict["external_magnetic"]["initializer"]}
        , b0xinit_{dict["external_magnetic"]["initializer"]["x_component"]
                       .template to<initializer::InitFunction<dimension>>()}
        , b0yinit_{dict["external_magnetic"]["initializer"]["y_component"]
                       .template to<initializer::InitFunction<dimension>>()}
        , b0zinit_{dict["external_magnetic"]["initializer"]["z_component"]
                       .template to<initializer::InitFunction<dimension>>()}
    {
    }

    //-------------------------------------------------------------------------
    //                  ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return B0.isUsable() and B0x_Ez.isUsable() and B0y_Ez.isUsable() and B0y_Ex.isUsable()
               and B0z_Ex.isUsable() and B0x_Ey.isUsable() and B0z_Ey.isUsable()
               and B0y_FaceX.isUsable() and B0z_FaceX.isUsable() and B0x_FaceY.isUsable()
               and B0z_FaceY.isUsable() and B0x_FaceZ.isUsable() and B0y_FaceZ.isUsable();
    }

    NO_DISCARD bool isSettable() const
    {
        return B0.isSettable() and B0x_Ez.isSettable() and B0y_Ez.isSettable()
               and B0y_Ex.isSettable() and B0z_Ex.isSettable() and B0x_Ey.isSettable()
               and B0z_Ey.isSettable() and B0y_FaceX.isSettable() and B0z_FaceX.isSettable()
               and B0x_FaceY.isSettable() and B0z_FaceY.isSettable() and B0x_FaceZ.isSettable()
               and B0y_FaceZ.isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(B0, B0x_Ez, B0y_Ez, B0y_Ex, B0z_Ex, B0x_Ey, B0z_Ey, B0y_FaceX,
                                     B0z_FaceX, B0x_FaceY, B0z_FaceY, B0x_FaceZ, B0y_FaceZ);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(B0, B0x_Ez, B0y_Ez, B0y_Ex, B0z_Ex, B0x_Ey, B0z_Ey, B0y_FaceX,
                                     B0z_FaceX, B0x_FaceY, B0z_FaceY, B0x_FaceZ, B0y_FaceZ);
    }

    //-------------------------------------------------------------------------

    // (Re-)evaluate B0 and all its analytic samples on a patch from the user functions. Called at
    // initialization and after every regrid; never interpolated across levels.
    template<typename GridLayout>
    void update(GridLayout const& layout)
    {
        B0init_.initialize(B0, layout);
        // EMF edges
        FieldUserFunctionInitializer::initialize(B0x_Ez, layout, b0xinit_);
        FieldUserFunctionInitializer::initialize(B0y_Ez, layout, b0yinit_);
        FieldUserFunctionInitializer::initialize(B0y_Ex, layout, b0yinit_);
        FieldUserFunctionInitializer::initialize(B0z_Ex, layout, b0zinit_);
        FieldUserFunctionInitializer::initialize(B0x_Ey, layout, b0xinit_);
        FieldUserFunctionInitializer::initialize(B0z_Ey, layout, b0zinit_);
        // Riemann transverse faces
        FieldUserFunctionInitializer::initialize(B0y_FaceX, layout, b0yinit_);
        FieldUserFunctionInitializer::initialize(B0z_FaceX, layout, b0zinit_);
        FieldUserFunctionInitializer::initialize(B0x_FaceY, layout, b0xinit_);
        FieldUserFunctionInitializer::initialize(B0z_FaceY, layout, b0zinit_);
        FieldUserFunctionInitializer::initialize(B0x_FaceZ, layout, b0xinit_);
        FieldUserFunctionInitializer::initialize(B0y_FaceZ, layout, b0yinit_);
    }

    // native (face-centered) background vector: B0x@Pdd, B0y@Dpd, B0z@Ddp
    VecFieldT B0;

    // analytic B0 samples at the constrained-transport EMF edges
    field_type B0x_Ez; // Ez edge (ppd)
    field_type B0y_Ez;
    field_type B0y_Ex; // Ex edge (dpp)
    field_type B0z_Ex;
    field_type B0x_Ey; // Ey edge (pdp)
    field_type B0z_Ey;

    // analytic B0 transverse samples at the Riemann face locations
    field_type B0y_FaceX; // centered like Bx (X-face)
    field_type B0z_FaceX;
    field_type B0x_FaceY; // centered like By (Y-face)
    field_type B0z_FaceY;
    field_type B0x_FaceZ; // centered like Bz (Z-face)
    field_type B0y_FaceZ;

private:
    VecFieldInitializer<dimension> B0init_;
    initializer::InitFunction<dimension> b0xinit_;
    initializer::InitFunction<dimension> b0yinit_;
    initializer::InitFunction<dimension> b0zinit_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_MAGNETIC_FIELD_HPP
