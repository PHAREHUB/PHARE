#ifndef PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/equality.hpp"
#include "core/hybrid/hybrid_quantities.hpp"

namespace PHARE::core
{

template<std::size_t dim>
using Field_t = Field<dim, HybridQuantity::Scalar, double>;


template<bool binary_eq = false>
struct FieldComparator
{
    auto float_eq(auto const a, auto const b) const
    {
        if constexpr (binary_eq)
            return a == b;
        else
            return any_float_eq(a, b, diff);
    };

    template<typename F0, typename F1>
    auto operator()(F0 const& ref, F1 const& cmp)
    {
        auto const& ref_dat = ref.data();
        auto const& cmp_dat = cmp.data();
        for (std::size_t i = 0; i < ref.size(); ++i)
        {
            ref0 = ref_dat[i] == 0 ? ref0 + 1 : ref0;
            cmp0 = cmp_dat[i] == 0 ? cmp0 + 1 : cmp0;

            nan0 = std::isnan(ref_dat[i]) ? nan0 + 1 : nan0;
            nan1 = std::isnan(cmp_dat[i]) ? nan1 + 1 : nan1;

            if (!(std::isnan(ref_dat[i]) || std::isnan(cmp_dat[i])))
            {
                auto ret = std::abs(ref_dat[i] - cmp_dat[i]);
                if (ret < diff)
                {
                    ++eqvals;
                    if (ref_dat[i] != 0 and cmp_dat[i] != 0)
                        ++eqnot0;
                }
                else
                    max_diff = ret > max_diff ? ret : max_diff;
            }
        }
        ok = eqvals == ref.size() and nan0 == 0 and nan1 == 0;
        return std::make_tuple(eqvals, eqnot0, ref0, cmp0);
    }

    operator bool() const { return ok; }

    double const diff  = 1e-15;
    std::size_t eqvals = 0, eqnot0 = 0, ref0 = 0, cmp0 = 0, nan0 = 0, nan1 = 0;
    bool ok         = true;
    double max_diff = 0;
};

using FloatFieldComparator_t = FieldComparator<false>;


template<std::size_t dim, typename PQ, typename D0, typename D1>
EqualityReport compare_fields(Field<dim, PQ, D0> const& ref, Field<dim, PQ, D1> const& cmp,
                              double const diff = 1e-15)
{
    auto const same_sizes = ref.size() == cmp.size();

    if (!same_sizes)
        return EqualityReport{false, "Tensorfield shape/size mismatch"};

    std::stringstream log;

    FloatFieldComparator_t eq{diff};
    auto const [eqvals, eqnot0, ref0, cmp0] = eq(ref, cmp);

    std::string const names
        = ref.name() == cmp.name() ? ref.name() : ref.name() + std::string{"/"} + cmp.name();
    log << "Fields compare (" << names << ") ";

    if (!eq)
    {
        auto const bad = ref.size() - eqvals;
        log << "value mismatch: \n";
        log << "ok(" << eqvals << ") - ";
        log << "ok!=0(" << eqnot0 << ") - ";
        log << "bad(" << bad << ") - ";
        log << "ref0(" << ref0 << ") - ";
        log << "cmp0(" << cmp0 << ") - ";
        log << "diff(" << eq.max_diff << ") - ";
        log << "nan0(" << eq.nan0 << ") - ";
        log << "nan1(" << eq.nan1 << ")\n";
        return EqualityReport{false, log.str()};
    }

    log << "are == with ";
    log << "ok(" << eqvals << ") - ";
    log << "ok!=0(" << eqnot0 << ")  ";

    return EqualityReport{true, log.str()};
}

template<typename... T0s, typename... T1s>
EqualityReport compare_fields(Grid<T0s...> const& ref, Grid<T1s...> const& cmp,
                              double const diff = 1e-15)
{
    return compare_fields(*ref, *cmp, diff);
}

} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP*/
