#ifndef DIAGNOSTIC_WRITER_H
#define DIAGNOSTIC_WRITER_H

#include "diagnostic_props.hpp"

#include "cppdict/include/dict.hpp"

namespace PHARE::diagnostic
{
class TypeWriter
{
public:
    virtual void write(DiagnosticProperties&)   = 0;
    virtual void compute(DiagnosticProperties&) = 0;
    virtual ~TypeWriter() {}
};


} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_WRITER_H
