#ifndef PHARE_CORE_UTILITIES_TIMESTAMPS_HPP
#define PHARE_CORE_UTILITIES_TIMESTAMPS_HPP

#include <string>
#include <cassert>
#include <cstdint>

#include "core/logger.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::core
{
struct ITimeStamper
{
    virtual double operator+=(double const& new_dt) noexcept = 0;

    virtual ~ITimeStamper() {}
};


class ConstantTimeStamper : public ITimeStamper
{
public:
    ConstantTimeStamper(double const& dt, std::size_t const& init_idx = 0)
        : dt_{dt}
        , idx_{init_idx}
    {
    }

    double operator+=([[maybe_unused]] double const& new_dt) noexcept override
    {
        assert(dt_ == new_dt); // binary comparison - should never fail in this case
        return dt_ * ++idx_;
    }

private:
    double dt_       = 0;
    std::size_t idx_ = 0;
};

struct TimeStamperFactory
{
    static std::unique_ptr<ITimeStamper> create(initializer::PHAREDict const& dict)
    {
        assert(dict.contains("time_step"));
        auto time_step  = dict["time_step"].template to<double>();
        std::size_t idx = 0;

        // only option for the moment https://github.com/PHAREHUB/PHARE/issues/475
        return std::make_unique<ConstantTimeStamper>(time_step, idx);
    }
};


} // namespace PHARE::core

#endif /*PHARE_CORE_UTILITIES_TIMESTAMPS_H */
