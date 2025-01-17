#ifndef PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP

#include <cassert>
#include <cmath>
#include <functional>

#include "phare_core.hpp"
#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

namespace PHARE::core
{

template<std::size_t dim>
class UsableElectromag : public Electromag<VecField_t<dim>>
{
    void _set()
    {
        E.set_on(Super::E);
        B.set_on(Super::B);
    }

public:
    using Super = Electromag<VecField_t<dim>>;

    template<typename GridLayout>
    UsableElectromag(GridLayout const& layout)
        : Super{"EM"}
        , E{"EM_E", layout, HybridQuantity::Vector::E}
        , B{"EM_B", layout, HybridQuantity::Vector::B}
    {
        _set();
    }

    UsableElectromag(UsableElectromag&& that)
        : Super{std::forward<Super>(that)}
        , E{std::move(that.E)}
        , B{std::move(that.B)}
    {
        _set();
    }

    Super& view() { return *this; }
    Super const& view() const { return *this; }
    auto& operator*() { return view(); }
    auto& operator*() const { return view(); }

    UsableVecField<dim> E, B;
};

} // namespace PHARE::core

#endif /* PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP */
