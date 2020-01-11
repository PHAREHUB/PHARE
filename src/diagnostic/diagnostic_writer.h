#ifndef DIAGNOSTIC_WRITER_H
#define DIAGNOSTIC_WRITER_H

#include "diagnostic_dao.h"

namespace PHARE::diagnostic
{
class DiagnosticTypeWriter
{
public:
    virtual void write(DiagnosticDAO&)   = 0;
    virtual void compute(DiagnosticDAO&) = 0;
    virtual ~DiagnosticTypeWriter() {}
};

class IDiagnosticWriter
{
public:
    virtual ~IDiagnosticWriter() {}
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_WRITER_H
