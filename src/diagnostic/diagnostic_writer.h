#ifndef DIAGNOSTIC_WRITER_H
#define DIAGNOSTIC_WRITER_H

#include "diagnostic_dao.h"

namespace PHARE::diagnostic
{
class DiagnosticWriter
{
public:
    virtual void write(DiagnosticDAO&)   = 0;
    virtual void compute(DiagnosticDAO&) = 0;
    virtual ~DiagnosticWriter() {}
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_WRITER_H
