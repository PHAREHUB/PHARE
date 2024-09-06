#ifndef _PHARE_AMR_DATA_ELECTROMAG_ELECTROMAG_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_ELECTROMAG_ELECTROMAG_INITIALIZER_HPP_

#include "core/data/vecfield/vecfield_initializer.hpp"
#include "amr/data/field/initializers/samrai_hdf5_field_initializer.hpp"
#include "initializer/data_provider.hpp"


#include <array>

namespace PHARE::amr
{

template<typename Electromag_t, typename GridLayout>
class ElectromagInitializer
{
public:
    virtual void init(Electromag_t& /*em*/, GridLayout const& /*layout*/) const { /*noop*/ }
    virtual ~ElectromagInitializer() {}
};

template<typename Electromag_t, typename GridLayout>
class ElectromagUserFuncInitializer : public ElectromagInitializer<Electromag_t, GridLayout>
{
    using VecFieldInit = core::VecFieldInitializer<Electromag_t::dimension>;

public:
    ElectromagUserFuncInitializer(initializer::PHAREDict const& dict)
        : Binit_{dict["magnetic"]["initializer"]}
    {
    }
    virtual ~ElectromagUserFuncInitializer() {}

    void init(Electromag_t& em, GridLayout const& layout) const override
    {
        Binit_.initialize(em.B, layout);
    }

    VecFieldInit Binit_;
};


/*
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/EM_B_y##default/d_box
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/EM_B_x##default/field_EM_B_x
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/EM_B_y##default/field_EM_B_y
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/EM_B_z##default/field_EM_B_z
*/


template<typename Electromag_t, typename GridLayout>
class ElectromagSamraiH5Initializer : public ElectromagInitializer<Electromag_t, GridLayout>
{
    using vecfield_type = typename Electromag_t::vecfield_type;
    using field_type    = typename vecfield_type::field_type;


public:
    ElectromagSamraiH5Initializer(initializer::PHAREDict const& dict)
        : dict_{dict}
    {
    }
    virtual ~ElectromagSamraiH5Initializer() {}

    void init(Electromag_t& em, GridLayout const& layout) const override
    {
        for (auto& field : em.B)
            SamraiHDF5FieldInitializer<field_type, GridLayout>{}.load(field, layout);
    }

    initializer::PHAREDict const dict_;
};



class ElectromagInitializerFactory
{
public:
    template<typename Electromag_t, typename GridLayout>
    NO_DISCARD static std::unique_ptr<ElectromagInitializer<Electromag_t, GridLayout>>
    create(initializer::PHAREDict const& dict)
    {
        if (dict["magnetic"]["initializer"].contains("x_component"))
            return std::make_unique<ElectromagUserFuncInitializer<Electromag_t, GridLayout>>(dict);
        else
            return std::make_unique<ElectromagSamraiH5Initializer<Electromag_t, GridLayout>>(dict);
    }
};

} // namespace PHARE::amr

#endif // _PHARE_AMR_DATA_ELECTROMAG_ELECTROMAG_INITIALIZER_HPP_
