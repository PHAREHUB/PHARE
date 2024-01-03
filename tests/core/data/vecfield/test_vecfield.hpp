

#ifndef PHARE_TEST_CORE_VECFIELD_TEST_HPP
#define PHARE_TEST_CORE_VECFIELD_TEST_HPP

#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield.hpp"

namespace PHARE::core
{
template<typename Field>
struct VecFieldMock
{
    using field_type                = Field;
    static auto constexpr dimension = Field::dimension;
    Field fm;
    Field& getComponent(Component) { return fm; }
    Field const& getComponent(Component) const { return fm; }

    auto& operator()(Component c) { return getComponent(c); }
    auto& operator()(Component c) const { return getComponent(c); }

    auto getComponents() const { return std::forward_as_tuple(fm, fm, fm); }
    auto getComponents() { return std::forward_as_tuple(fm, fm, fm); }

    auto operator()() { return getComponents(); }
    auto operator()() const { return getComponents(); }

    bool isUsable() const { return true; }


    auto& view() { return *this; }
    auto view() const { return *this; }
};


} // namespace PHARE::core



#endif /*PHARE_TEST_CORE_VECFIELD_TEST_HPP*/
