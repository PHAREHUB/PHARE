#ifndef DIAGNOSTIC_WRITER_HPP
#define DIAGNOSTIC_WRITER_HPP

#include "diagnostic_props.hpp"



namespace PHARE::diagnostic
{
class TypeWriter
{
public:
    virtual void write(DiagnosticProperties&)   = 0;
    virtual void compute(DiagnosticProperties&) = 0;

    // request the model temporaries the diagnostic needs, when the diagnostic is declared
    virtual void requestTemporaries(DiagnosticProperties&) {}

    virtual ~TypeWriter() {}
};


} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_WRITER_HPP
