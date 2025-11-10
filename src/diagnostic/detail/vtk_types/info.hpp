#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_INFO_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_INFO_HPP

#include "diagnostic/detail/vtkh5_type_writer.hpp"


namespace PHARE::diagnostic::vtkh5
{


template<typename H5Writer>
class InfoDiagnosticWriter : public H5TypeWriter<H5Writer>
{
    using Super         = H5TypeWriter<H5Writer>;
    using VTKFileWriter = Super::VTKFileWriter;

public:
    InfoDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void write(DiagnosticProperties&) override;

    void compute(DiagnosticProperties&) override {}
};



template<typename H5Writer>
void InfoDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
}


} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_INFO_H */
