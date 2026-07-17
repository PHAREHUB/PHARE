#ifndef PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_HPP
#define PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_HPP

#include <memory>
#include <vector>

#include "core/utilities/span.hpp"

namespace PHARE::initializer::test_fn::func_1d
{
using Param  = std::vector<double> const&;
using Return = std::shared_ptr<PHARE::core::Span<double>>;

Return inline density(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x.size(), 1);
}

Return inline vx(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vy(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vz(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthx(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthy(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthz(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline bx(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline by(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline bz(Param x)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}


} // namespace PHARE::initializer::test_fn::func_1d


namespace PHARE::initializer::test_fn::func_2d
{
using Param  = std::vector<double> const&;
using Return = std::shared_ptr<PHARE::core::Span<double>>;

Return inline density(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x.size(), 1);
}

Return inline vx(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vy(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vz(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthx(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthy(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthz(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline bx(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline by(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline bz(Param x, Param /*y*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}


} // namespace PHARE::initializer::test_fn::func_2d

namespace PHARE::initializer::test_fn::func_3d
{
using Param  = std::vector<double> const&;
using Return = std::shared_ptr<PHARE::core::Span<double>>;

Return inline density(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x.size(), 1);
}

Return inline vx(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vy(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vz(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthx(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthy(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline vthz(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline bx(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline by(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return inline bz(Param x, Param /*y*/, Param /*z*/)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}


} // namespace PHARE::initializer::test_fn::func_3d

template<std::size_t dim>
auto makeSharedPtr()
{
    using Param = std::vector<double> const&;

    if constexpr (dim == 1)
    {
        return [](Param x) { return std::make_shared<PHARE::core::VectorSpan<double>>(x); };
    }
    else if constexpr (dim == 2)
    {
        return [](Param x, Param /*y*/) {
            return std::make_shared<PHARE::core::VectorSpan<double>>(x);
        };
    }
    else if constexpr (dim == 3)
    {
        return [](Param x, Param /*y*/, Param /*z*/) {
            return std::make_shared<PHARE::core::VectorSpan<double>>(x);
        };
    }
}


#endif // PHARE_TEST_INITIALIZER_INIT_FUNCTIONS_HPP
