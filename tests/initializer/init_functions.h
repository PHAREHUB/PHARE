#ifndef PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_H
#define PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_H

#include <memory>
#include <vector>

#include "core/utilities/span.h"

namespace PHARE::initializer::test_fn::func_1d
{
using Param  = std::vector<double> const&;
using Return = std::shared_ptr<PHARE::core::Span<double>>;

Return density(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vx(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vy(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vz(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vthx(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vthy(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vthz(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return bx(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return by(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return bz(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

} // namespace PHARE::initializer::test_fn::func_1d


namespace PHARE::initializer::test_fn::func_2d
{
using Param  = std::vector<double> const&;
using Return = std::shared_ptr<PHARE::core::Span<double>>;

Return density(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vx(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vy(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vz(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vthx(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vthy(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vthz(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return bx(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return by(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return bz(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

} // namespace PHARE::initializer::test_fn::func_2d

#endif // PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_H
