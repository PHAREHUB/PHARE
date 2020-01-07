#ifndef PHARE_TEST_DIAGNOSTIC_INCLUDE_H
#define PHARE_TEST_DIAGNOSTIC_INCLUDE_H

#include "tests/simulator/per_test.h"

#include "diagnostic/detail/highfive.h"
#include "diagnostic/detail/types/electromag.h"
#include "diagnostic/detail/types/particle.h"
#include "diagnostic/detail/types/fluid.h"

using namespace PHARE::diagnostic;
using namespace PHARE::diagnostic::h5;

constexpr unsigned NEW_HI5_FILE
    = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate;

template<typename Hierarchy, typename HybridModel>
struct Hi5Diagnostic
{
    // using HybridModelT        = typename Simulator::HybridModel;
    using DiagnosticModelView = AMRDiagnosticModelView<Hierarchy, HybridModel>;
    using DiagnosticWriter    = HighFiveDiagnostic<DiagnosticModelView>;

    Hi5Diagnostic(Hierarchy& hierarchy, HybridModel& hybridModel, std::string fileName,
                  unsigned flags)
        : hierarchy_{hierarchy}
        , model_{hybridModel}
        , file_{fileName}
        , flags_{flags}
    {
    }
    ~Hi5Diagnostic() {}

    auto dict(std::string&& type, std::string& subtype)
    {
        PHARE::initializer::PHAREDict dict;
        dict["name"]            = type;
        dict["type"]            = type;
        dict["subtype"]         = subtype;
        dict["compute_every"]   = std::size_t{1};
        dict["write_every"]     = std::size_t{1};
        dict["start_iteration"] = std::size_t{0};
        dict["last_iteration"]  = std::numeric_limits<std::size_t>::max();
        return dict;
    }
    auto electromag(std::string&& subtype) { return dict("electromag", subtype); }
    auto particles(std::string&& subtype) { return dict("particles", subtype); }
    auto fluid(std::string&& subtype) { return dict("fluid", subtype); }

    std::string filename() const { return file_ + "_hi5_test.5"; }

    std::string getPatchPath(int level, PHARE::amr::SAMRAI_Types::patch_t& patch)
    {
        std::stringstream globalId;
        globalId << patch.getGlobalId();
        return writer.getPatchPath("time", level, globalId.str());
    }

    template<typename Field>
    void checkField(Field& field, std::string path)
    {
        std::vector<float> fieldV;
        writer.file().getDataSet(path).read(fieldV);
        EXPECT_EQ(fieldV.size(), field.size());

        for (size_t i = 0; i < fieldV.size(); i++)
        {
            if (!std::isnan(fieldV[i]))
            {
                EXPECT_FLOAT_EQ(fieldV[i], field.data()[i]);
            }
        }
    }

    template<typename VecField>
    void checkVecField(VecField& vecField, std::string fieldPath)
    {
        for (auto& [id, type] : PHARE::core::Components::componentMap)
        {
            checkField(vecField.getComponent(type), fieldPath + "/" + id);
        }
    }


    Hierarchy& hierarchy_;
    HybridModel& model_;
    std::string file_;
    unsigned flags_;
    DiagnosticModelView modelView{hierarchy_, model_};
    DiagnosticWriter writer{modelView, filename(), flags_};
    DiagnosticsManager<DiagnosticWriter> dMan{writer};
};


#endif /*PHARE_TEST_DIAGNOSTIC_INCLUDE_H*/
