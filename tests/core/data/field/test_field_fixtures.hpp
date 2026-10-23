#ifndef PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/equality.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"


namespace PHARE::core
{

template<typename core_types>
struct TestFieldOptions
{
    auto constexpr static opts = core_types::opts;

    using GridLayout_t = core_types::GridLayout_t;
};


// the outermost ghost node of some quantities never receives a valid
// interpolated value (last ghost can have no interpolation), so it must be
// excluded from comparisons that otherwise cover the ghost box.
template<typename PQ>
bool valid_ghost_box(PQ const physicalQuantity)
{
    if constexpr (std::is_same_v<PQ, HybridQuantity::Scalar>)
    {
        using enum HybridQuantity::Scalar;
        return not any_in(physicalQuantity, rho, Vx, Vy, Vz);
    }
    throw std::runtime_error("No other impl");
}


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

    void add(auto const a, auto const b)
    {
        ++total;
        ref0 = a == 0 ? ref0 + 1 : ref0;
        cmp0 = b == 0 ? cmp0 + 1 : cmp0;

        nan0 = std::isnan(a) ? nan0 + 1 : nan0;
        nan1 = std::isnan(b) ? nan1 + 1 : nan1;

        if (!(std::isnan(a) || std::isnan(b)))
        {
            auto ret = std::abs(a - b);
            if (ret < diff)
            {
                ++eqvals;
                if (a != 0 and b != 0)
                    ++eqnot0;
            }
            else
                max_diff = ret > max_diff ? ret : max_diff;
        }
    }

    void finalize() { ok = eqvals == total and nan0 == 0 and nan1 == 0; }

    template<typename F0, typename F1>
    auto operator()(F0 const& ref, F1 const& cmp)
    {
        auto const& ref_dat = ref.data();
        auto const& cmp_dat = cmp.data();
        for (std::size_t i = 0; i < ref.size(); ++i)
            add(ref_dat[i], cmp_dat[i]);
        finalize();
        return std::make_tuple(eqvals, eqnot0, ref0, cmp0);
    }

    operator bool() const { return ok; }

    double const diff  = 1e-15;
    std::size_t total  = 0;
    std::size_t eqvals = 0, eqnot0 = 0, ref0 = 0, cmp0 = 0, nan0 = 0, nan1 = 0;
    bool ok         = true;
    double max_diff = 0;
};

using FloatFieldComparator_t = FieldComparator<false>;


template<std::size_t dim, typename PQ, typename D0, typename D1, auto am0, auto am1>
EqualityReport compare_fields(Field<dim, PQ, D0, am0> const& ref,
                              Field<dim, PQ, D1, am1> const& cmp, double const diff = 1e-15)
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


template<typename... T0s, typename... T1s>
EqualityReport compare_fields(FieldTileSet<T0s...> const& ref, FieldTileSet<T1s...> const& cmp,
                              [[maybe_unused]] double const diff = 1e-15)
{
    auto const same_sizes = ref.size() == cmp.size() and ref().size() == cmp().size();
    if (!same_sizes)
        return EqualityReport{false, "Tensorfield shape/size mismatch"};

    std::string const names
        = ref.name() == cmp.name() ? ref.name() : ref.name() + std::string{"/"} + cmp.name();
    std::stringstream log;
    log << "Fields compare (" << names << ") ";

    std::size_t ok = 0, okn0 = 0;
    for (std::size_t tidx = 0; tidx < ref().size(); ++tidx)
    {
        FloatFieldComparator_t eq{diff};
        auto const [eqvals, eqnot0, ref0, cmp0] = eq(ref()[tidx](), cmp()[tidx]());

        if (!eq)
        {
            auto const bad = ref.size() - eqvals;
            log << "value mismatch: \n";
            log << "ok(" << eqvals << ") - ";
            log << "ok!=0(" << eqnot0 << ") - ";
            log << "bad(" << bad << ") - ";
            log << "ref0(" << ref0 << ") - ";
            log << "cmp0(" << cmp0 << ") - ";
            log << "nan0(" << eq.nan0 << ") - ";
            log << "nan1(" << eq.nan1 << ")\n";
            return EqualityReport{false, log.str()};
        }

        ok += eqvals;
        okn0 += eqnot0;
    }

    log << "are == with ";
    log << "ok(" << ok << ") - ";
    log << "ok!=0(" << okn0 << ")  ";

    return EqualityReport{true, log.str()};
}

template<typename... T0s, typename... T1s>
EqualityReport compare_fields(GridTileSet<T0s...> const& ref, GridTileSet<T1s...> const& cmp,
                              [[maybe_unused]] double const diff = 1e-15)
{
    return compare_field(*ref, *cmp, diff);
}


template<std::size_t dim, typename PQ, typename D, auto am, typename... T1s>
EqualityReport compare_fields(Field<dim, PQ, D, am> const& ref, FieldTileSet<T1s...> const& cmp,
                              double const diff = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "Field/FieldTileSet shape mismatch"};
    auto tmp = ref;
    tmp.zero();
    reduce_single(tmp, cmp);
    return compare_fields(ref, tmp, diff);
}

template<typename... T0s, std::size_t dim, typename PQ, typename D, auto am>
EqualityReport compare_fields(FieldTileSet<T0s...> const& ref, Field<dim, PQ, D, am> const& cmp,
                              double const diff = 1e-15)
{
    return compare_fields(cmp, ref, diff);
}


/**
 * @brief compare_field_domains is compare_fields() restricted to the domain box, for use
 * on quantities where the ghost box is not valid_ghost_box() (i.e. some ghost nodes are
 * never given a valid interpolated value, e.g. density and bulk velocity).
 * Supports both AoSMapped (plain Field) and AoSPCTS (FieldTileSet) inputs.
 */
template<typename GridLayout, std::size_t dim, typename PQ, typename D0, typename D1, auto am0,
        auto am1>
EqualityReport compare_field_domains(GridLayout const& layout, Field<dim, PQ, D0, am0> const& ref,
                                     Field<dim, PQ, D1, am1> const& cmp,
                                     double const diff = 1e-15)
{
    auto const same_sizes = ref.size() == cmp.size();

    if (!same_sizes)
        return EqualityReport{false, "Tensorfield shape/size mismatch"};

    FloatFieldComparator_t eq{diff};
    layout.evalOnBox(ref, [&](auto const&... idxs) { eq.add(ref(idxs...), cmp(idxs...)); });
    eq.finalize();

    std::stringstream log;
    std::string const names
        = ref.name() == cmp.name() ? ref.name() : ref.name() + std::string{"/"} + cmp.name();
    log << "Fields compare (domain) (" << names << ") ";

    if (!eq)
    {
        auto const bad = eq.total - eq.eqvals;
        log << "value mismatch: \n";
        log << "ok(" << eq.eqvals << ") - ";
        log << "ok!=0(" << eq.eqnot0 << ") - ";
        log << "bad(" << bad << ") - ";
        log << "ref0(" << eq.ref0 << ") - ";
        log << "cmp0(" << eq.cmp0 << ") - ";
        log << "diff(" << eq.max_diff << ") - ";
        log << "nan0(" << eq.nan0 << ") - ";
        log << "nan1(" << eq.nan1 << ")\n";
        return EqualityReport{false, log.str()};
    }

    log << "are == with ";
    log << "ok(" << eq.eqvals << ") - ";
    log << "ok!=0(" << eq.eqnot0 << ")  ";

    return EqualityReport{true, log.str()};
}

template<typename GridLayout, typename... T0s, typename... T1s>
EqualityReport compare_field_domains(GridLayout const& layout, Grid<T0s...> const& ref,
                                     Grid<T1s...> const& cmp, double const diff = 1e-15)
{
    return compare_field_domains(layout, *ref, *cmp, diff);
}

template<typename GridLayout, std::size_t dim, typename PQ, typename D, auto am, typename... T1s>
EqualityReport compare_field_domains(GridLayout const& layout, Field<dim, PQ, D, am> const& ref,
                                     FieldTileSet<T1s...> const& cmp, double const diff = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "Field/FieldTileSet shape mismatch"};
    auto tmp = ref;
    tmp.zero();
    reduce_single(tmp, cmp);
    return compare_field_domains(layout, ref, tmp, diff);
}

template<typename GridLayout, typename... T0s, std::size_t dim, typename PQ, typename D, auto am>
EqualityReport compare_field_domains(GridLayout const& layout, FieldTileSet<T0s...> const& ref,
                                     Field<dim, PQ, D, am> const& cmp, double const diff = 1e-15)
{
    return compare_field_domains(layout, cmp, ref, diff);
}


} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_HPP*/
