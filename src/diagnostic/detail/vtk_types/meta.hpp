#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_META_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_META_HPP

#include "diagnostic/detail/vtkh5_type_writer.hpp"


namespace PHARE::diagnostic::vtkh5
{


template<typename H5Writer>
class MetaDiagnosticWriter : public H5TypeWriter<H5Writer>
{
    using Super         = H5TypeWriter<H5Writer>;
    using VTKFileWriter = Super::VTKFileWriter;

public:
    MetaDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void write(DiagnosticProperties&) override;

    void compute(DiagnosticProperties&) override {}
};

template<typename H5Writer>
void MetaDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
}



} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_META_HPP */
