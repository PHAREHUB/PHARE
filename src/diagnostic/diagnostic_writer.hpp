#ifndef DIAGNOSTIC_WRITER_HPP
#define DIAGNOSTIC_WRITER_HPP

#include "diagnostic_props.hpp"



namespace PHARE::diagnostic
{
class TypeWriter
{
public:
    virtual void write(DiagnosticProperties&) = 0;

    // called by the writer immediately before write(), computes into the model fields
    virtual void compute_as_needed(DiagnosticProperties&) = 0;

    // called by the DiagnosticsManager on compute_timestamps. Unused today: no writer overrides
    // it and pharein sets compute_timestamps to write_timestamps. Kept for computations that
    // would run more often than writes (e.g. time averaging)
    virtual void compute(DiagnosticProperties&) {}

    virtual ~TypeWriter() {}
};


} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_WRITER_HPP
