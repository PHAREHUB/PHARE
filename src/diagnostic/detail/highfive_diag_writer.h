#ifndef HIGHFIVEDIAGNOSTICWRITER_H
#define HIGHFIVEDIAGNOSTICWRITER_H

#include <string>


#include "diagnostic_writer.h"

namespace PHARE
{
namespace diagnostic
{
    namespace h5
    {
        template<typename HighFiveDiagnostic>
        class Hi5DiagnosticWriter : public DiagnosticWriter
        {
        public:
            using Attributes = typename HighFiveDiagnostic::Attributes;
            Hi5DiagnosticWriter(HighFiveDiagnostic& hi5)
                : hi5_(hi5)
            {
            }

            virtual void getDataSetInfo(std::string const& patchID, Attributes& patchAttributes)
                = 0;
            virtual void initDataSets(std::vector<std::string> const& patchIDs,
                                      Attributes& patchAttributes){};

        protected:
            HighFiveDiagnostic& hi5_;
        };
    } // namespace h5
} // namespace diagnostic
} // namespace PHARE


#endif // HIGHFIVEDIAGNOSTICWRITER_H
