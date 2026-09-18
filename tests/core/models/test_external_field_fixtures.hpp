#ifndef PHARE_TEST_CORE_MODELS_TEST_EXTERNAL_FIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_MODELS_TEST_EXTERNAL_FIELD_FIXTURES_HPP

#include "core/models/external_field.hpp"
#include "core/utilities/point/point.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/span.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace PHARE::core
{
/**
 * @brief turn a point-wise formula f(Point, time) into a vectorized SpaceTimeFunction.
 *
 * Stands in for pyphare's space_time_fn_wrapper: coordinates arrive as Span views
 * and the result is returned as a Span owning its buffer, so a test drives the updaters
 * through the very call convention the python binding produces, without python.
 */
template<std::size_t dim, typename Fn>
initializer::SpaceTimeFunction<dim> spaceTimeFunction(Fn f)
{
    auto fill = [f](Span<double const> const& x, double t, auto&& pointAt) {
        std::vector<double> out(x.size());
        for (std::size_t i = 0; i < x.size(); ++i)
            out[i] = f(pointAt(i), t);
        return std::static_pointer_cast<Span<double>>(
            std::make_shared<VectorSpan<double>>(std::move(out)));
    };

    if constexpr (dim == 1)
        return [fill](Span<double const> const& x, double t) {
            return fill(x, t, [&](std::size_t i) { return Point<double, 1>{x[i]}; });
        };
    else if constexpr (dim == 2)
        return [fill](Span<double const> const& x, Span<double const> const& y, double t) {
            return fill(x, t, [&](std::size_t i) { return Point<double, 2>{x[i], y[i]}; });
        };
    else
        return [fill](Span<double const> const& x, Span<double const> const& y,
                      Span<double const> const& z, double t) {
            return fill(x, t, [&](std::size_t i) { return Point<double, 3>{x[i], y[i], z[i]}; });
        };
}


/**
 * @brief an ExternalField that owns the memory of its vecfields, for tests.
 *
 * Mirrors what the resources manager does on a patch: the fields of an ExternalField are views,
 * whose buffers are set here from grids owned by this fixture.
 */
template<std::size_t dim>
class UsableExternalField : public ExternalField<VecFieldMHD<dim>>
{
public:
    using Super = ExternalField<VecFieldMHD<dim>>;

    template<typename GridLayout>
    UsableExternalField(std::string const& name, GridLayout const& layout)
        : Super{name}
        , b0_{name + "_B0", layout, MHDQuantity::Vector::B}
        , dB0dt_{name + "_dB0dt", layout, MHDQuantity::Vector::B}
        , tmpVec_{name + "_tmpVec", layout, MHDQuantity::Vector::VecAllPrimal}
        , scratch_{view_as(tmpVec_.super(), MHDQuantity::Vector::E, layout)}
    {
        b0_.set_on(this->B0);
        dB0dt_.set_on(this->dB0dt);
    }

    Super& super() { return *this; }
    VecFieldMHD<dim>& scratch() { return scratch_; }

private:
    UsableVecFieldMHD<dim> b0_;
    UsableVecFieldMHD<dim> dB0dt_;
    UsableVecFieldMHD<dim> tmpVec_;
    VecFieldMHD<dim> scratch_;
};

} // namespace PHARE::core

#endif /* PHARE_TEST_CORE_MODELS_TEST_EXTERNAL_FIELD_FIXTURES_HPP */
