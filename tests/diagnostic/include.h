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

template<typename GridLayout, typename Field, typename FieldFilter = PHARE::FieldNullFilter>
void checkField(HighFive::File& file, GridLayout& layout, Field& field, std::string path,
                FieldFilter ff = FieldFilter{})
{
    std::vector<float> fieldV;
    file.getDataSet(path).read(fieldV);
    EXPECT_EQ(fieldV.size(), field.size());
    size_t start = ff.template start<Field>(layout, PHARE::core::Direction::X);
    size_t end   = ff.template end<Field>(layout, PHARE::core::Direction::X);
    for (size_t i = start; i < end; i++)
    {
        if (std::isnan(fieldV[i]) || std::isnan(field.data()[i]))
            throw std::runtime_error("This field should not be NaN");

        EXPECT_FLOAT_EQ(fieldV[i], field.data()[i]);
    }
}


template<typename GridLayout, typename VecField, typename FieldFilter = PHARE::FieldNullFilter>
void checkVecField(HighFive::File& file, GridLayout& layout, VecField& vecField, std::string path,
                   FieldFilter ff = FieldFilter{})
{
    for (auto& [id, type] : PHARE::core::Components::componentMap)
    {
        checkField(file, layout, vecField.getComponent(type), path + "/" + id, ff);
    }
}


template<typename Hierarchy, typename HybridModel>
struct Hi5Diagnostic
{
    using DiagnosticModelView = AMRDiagnosticModelView<Hierarchy, HybridModel>;
    using DiagnosticWriter    = HighFiveDiagnosticWriter<DiagnosticModelView>;

    Hi5Diagnostic(Hierarchy& hierarchy, HybridModel& hybridModel, unsigned flags)
        : hierarchy_{hierarchy}
        , model_{hybridModel}
        , flags_{flags}
    {
    }
    ~Hi5Diagnostic() {}

    auto dict(std::string&& category, std::string& type)
    {
        PHARE::initializer::PHAREDict dict;
        dict["name"]            = type;
        dict["category"]        = category;
        dict["type"]            = type;
        dict["compute_every"]   = std::size_t{1};
        dict["write_every"]     = std::size_t{1};
        dict["start_iteration"] = std::size_t{0};
        dict["last_iteration"]  = std::numeric_limits<std::size_t>::max();
        return dict;
    }
    auto electromag(std::string&& type) { return dict("electromag", type); }
    auto particles(std::string&& type) { return dict("particle", type); }
    auto fluid(std::string&& type) { return dict("fluid", type); }

    std::string getPatchPath(int level, PHARE::amr::SAMRAI_Types::patch_t& patch)
    {
        std::stringstream globalId;
        globalId << patch.getGlobalId();
        return writer.getPatchPath("time", level, globalId.str());
    }



    Hierarchy& hierarchy_;
    HybridModel& model_;
    unsigned flags_;
    DiagnosticModelView modelView{hierarchy_, model_};
    DiagnosticWriter writer{modelView, "phare_outputs", flags_};
    DiagnosticsManager<DiagnosticWriter> dMan{writer};
};


#endif /*PHARE_TEST_DIAGNOSTIC_INCLUDE_H*/
