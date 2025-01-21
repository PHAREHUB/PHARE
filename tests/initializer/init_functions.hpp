#ifndef PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_HPP
#define PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_HPP

#include <memory>
#include <vector>

#include "core/utilities/span.hpp"


namespace PHARE::initializer::test_fn::func_1d
{

using namespace PHARE::core;

using Param  = std::vector<floater_t<4>> const&;
using Return = std::shared_ptr<PHARE::core::Span<floater_t<4>>>;

Return density(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vx(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vy(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vz(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vthx(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vthy(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vthz(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return bx(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return by(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return bz(Param x)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}


} // namespace PHARE::initializer::test_fn::func_1d


namespace PHARE::initializer::test_fn::func_2d
{

using namespace PHARE::core;

using Param  = std::vector<floater_t<4>> const&;
using Return = std::shared_ptr<PHARE::core::Span<floater_t<4>>>;

Return density(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vx(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vy(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vz(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vthx(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vthy(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return vthz(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return bx(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return by(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}

Return bz(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<floater_t<4>>>(x);
}


} // namespace PHARE::initializer::test_fn::func_2d


template<std::size_t dim>
auto makeSharedPtr()
{
    using namespace PHARE::core;

    using Param = std::vector<floater_t<4>> const&;

    if constexpr (dim == 1)
    {
        return [](Param x) { return std::make_shared<PHARE::core::VectorSpan<floater_t<4>>>(x); };
    }
    else if constexpr (dim == 2)
    {
        return [](Param x, Param /*y*/) {
            return std::make_shared<PHARE::core::VectorSpan<floater_t<4>>>(x);
        };
    }
    else if constexpr (dim == 3)
    {
        return [](Param x, Param /*y*/, Param /*z*/) {
            return std::make_shared<PHARE::core::VectorSpan<floater_t<4>>>(x);
        };
    }
}


#endif // PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_HPP
