#ifndef DIAGNOSTIC_WRITER_H
#define DIAGNOSTIC_WRITER_H

#include "diagnostic_props.h"

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

class IWriter
{
public:
    virtual ~IWriter() {}
    virtual void dump(std::vector<DiagnosticProperties*> const& diagnostics, double timestamp) = 0;
    virtual void dump_level(std::size_t level,
                            std::vector<DiagnosticProperties*> const& diagnostics, double timestamp)
        = 0;
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_WRITER_H
