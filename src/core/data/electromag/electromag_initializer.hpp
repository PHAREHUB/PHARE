#ifndef _PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_INITIALIZER_HPP_
#define _PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_INITIALIZER_HPP_

#include "initializer/data_provider.hpp"

#include <array>


namespace PHARE::core
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

    VecFieldInitializer<Electromag_t::dimension> Binit_;
};


template<typename Electromag_t, typename GridLayout>
class ElectromagInitializerFactory
{
public:
    NO_DISCARD static std::unique_ptr<ElectromagInitializer<Electromag_t, GridLayout>>
    create(initializer::PHAREDict const& dict)
    {
        return std::make_unique<ElectromagUserFuncInitializer<Electromag_t, GridLayout>>(dict);
        // else
        // return std::make_unique<ElectromagInitializer<Electromag_t>>();
    }
};

} // namespace PHARE::core

#endif // _PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_INITIALIZER_HPP_
