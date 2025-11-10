#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_PARTICLE_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_PARTICLE_HPP

#include "diagnostic/detail/vtkh5_type_writer.hpp"


// not sure this makes sense

namespace PHARE::diagnostic::vtkh5
{

template<typename H5Writer>
class ParticlesDiagnosticWriter : public H5TypeWriter<H5Writer>
{
    using Super         = H5TypeWriter<H5Writer>;
    using VTKFileWriter = Super::VTKFileWriter;

public:
    ParticlesDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {}
};



template<typename H5Writer>
void ParticlesDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
}



} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H */
