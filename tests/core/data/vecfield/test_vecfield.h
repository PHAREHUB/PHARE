

#ifndef PHARE_TEST_CORE_VECFIELD_TEST_H
#define PHARE_TEST_CORE_VECFIELD_TEST_H

#include "core/data/field/field.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield.h"

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


    auto& as_view() { return *this; }
    auto as_view() const { return *this; }
};


} // namespace PHARE::core



#endif /*PHARE_TEST_CORE_VECFIELD_TEST_H*/
