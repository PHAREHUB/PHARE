#ifndef PHARE_DIAGNOSTIC_HI5_DIAGNOSTIC_HPP_
#define PHARE_DIAGNOSTIC_HI5_DIAGNOSTIC_HPP_

#if !defined(PHARE_WITH_HIGHFIVE)
#define PHARE_WITH_HIGHFIVE
#endif
#include "phare/diagnostic_manager.h"

namespace PHARE
{
namespace hi5
{
    struct Diagnostic
    {
        HighFive::File& file_;
        PHARE::diagnostic::Mode mode_ = PHARE::diagnostic::Mode::LIGHT;

        Diagnostic(const Diagnostic&)             = delete;
        Diagnostic(const Diagnostic&&)            = delete;
        Diagnostic& operator&(const Diagnostic&)  = delete;
        Diagnostic& operator&(const Diagnostic&&) = delete;
    };

} // namespace hi5
} // namespace PHARE

#endif //  PHARE_DIAGNOSTIC_HI5_DIAGNOSTIC_HPP_