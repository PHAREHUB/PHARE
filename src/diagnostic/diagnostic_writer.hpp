#ifndef DIAGNOSTIC_WRITER_HPP
#define DIAGNOSTIC_WRITER_HPP

#include "diagnostic_props.hpp"



namespace PHARE::diagnostic
{
class TypeWriter
{
public:
    virtual void write(DiagnosticProperties&) = 0;

    // called by the writer immediately before write(), computes to temporaries
    virtual void compute_as_needed(DiagnosticProperties&) = 0;

    // called by the DiagnosticsManager on compute_timestamps, which may be more frequent than
    // write_timestamps (e.g. for time averaging)
    virtual void compute(DiagnosticProperties&) {}

    virtual ~TypeWriter() {}
};


} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_WRITER_HPP
