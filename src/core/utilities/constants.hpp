#ifndef PHARE_CORE_UTILITIES_CONSTANTS_HPP
#define PHARE_CORE_UTILITIES_CONSTANTS_HPP

#include <cstdint>
#include <cstddef>
#include <type_traits>
namespace PHARE
{
namespace core
{
    constexpr std::uint32_t NBR_COMPO{3};

    constexpr std::size_t dirX{0};
    constexpr std::size_t dirY{1};
    constexpr std::size_t dirZ{2};

    template<std::size_t Constant>
    class StrongIntegralConstant
    {
    public:
        constexpr decltype(auto) operator()() const { return constant(); }

    protected:
        static constexpr std::integral_constant<std::size_t, Constant> constant{};
    };

    template<std::size_t Constant>
    class DimConst : public StrongIntegralConstant<Constant>
    {
    };
    template<std::size_t Constant>
    class InterpConst : public StrongIntegralConstant<Constant>
    {
    };
    template<std::size_t Constant>
    class RefinedParticlesConst : public StrongIntegralConstant<Constant>
    {
    };

} // namespace core
} // namespace PHARE

#endif
