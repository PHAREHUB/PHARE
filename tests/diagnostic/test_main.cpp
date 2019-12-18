

#include <mpi.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "include.h"
#include "samrai_lifecycle.h"

#include "func.h"
#include "defs.h"
using namespace PHARE_test::_1d;

#include "integrator.h"
#include "tag_strat.h"
#include "hierarchy.h"

#include "diagnostic/detail/highfive.h"
#include "diagnostic/detail/types/electromag.h"
#include "diagnostic/detail/types/particle.h"
#include "diagnostic/detail/types/fluid.h"

using namespace PHARE::diagnostic::h5;
using namespace PHARE::diagnostic;


struct Hi5Diagnostic : public ::testing::Test, public AfullHybridBasicHierarchy
{
    ~Hi5Diagnostic() {}
    Hi5Diagnostic(std::string fileName,
      unsigned flags = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate) : file_{fileName}, flags_{flags}{}

    auto dict(std::string&& type, std::string& subtype)
    {
        PHARE::initializer::PHAREDict dict;
        dict["name"]            = type;
        dict["type"]            = type;
        dict["subtype"]         = subtype;
        dict["compute_every"]   = std::size_t{1};
        dict["write_every"]     = std::size_t{1};
        dict["start_iteration"] = std::size_t{0};
        dict["end_iteration"]   = std::numeric_limits<std::size_t>::max();
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
        for (auto& [id, type] : Components::componentMap)
        {
            checkField(vecField.getComponent(type), fieldPath + "/" + id);
        }
    }

    using DiagnosticModelView = AMRDiagnosticModelView<PHARE::amr::SAMRAI_Types, HybridModelT>;
    using DiagnosticWriter    = HighFiveDiagnostic<DiagnosticModelView>;

    std::string file_;
    unsigned flags_;
    DiagnosticModelView modelView{basicHierarchy->getHierarchy(), *hybridModel};
    DiagnosticWriter writer{modelView, filename(), flags_};
    DiagnosticsManager<DiagnosticWriter> dMan{writer};
};

struct FluidHi5Diagnostic : public Hi5Diagnostic{
    FluidHi5Diagnostic() : Hi5Diagnostic{"fluid"}{}
};
TEST_F(FluidHi5Diagnostic, fluidWrite)
{
    dMan.addDiagDict(this->fluid("/ions/density"))
        .addDiagDict(this->fluid("/ions/bulkVelocity"))
        .addDiagDict(this->fluid("/ions/pop/ions_alpha/density"))
        .addDiagDict(this->fluid("/ions/pop/ions_alpha/flux"))
        .addDiagDict(this->fluid("/ions/pop/ions_protons/density"))
        .addDiagDict(this->fluid("/ions/pop/ions_protons/flux"))
        .dump();
}

struct ReadOnlyFluidHi5Diagnostic : public Hi5Diagnostic{
    ReadOnlyFluidHi5Diagnostic() : Hi5Diagnostic{"fluid", HighFive::File::ReadOnly}{}
};
TEST_F(ReadOnlyFluidHi5Diagnostic, fluidRead)
{
    auto checkIons = [&](auto& ions, auto patchPath) {
        std::string path(patchPath + "/ions/");

        for (auto& pop : modelView.getIons())
        {
            std::string popPath(path + "pop/" + pop.name() + "/");
            checkField(pop.density(), popPath + "density");
            checkVecField(pop.flux(), popPath + "flux");
        }

        checkField(ions.density(), path + "density");
        checkVecField(ions.velocity(), path + "bulkVelocity");
    };

    auto& hierarchy = basicHierarchy->getHierarchy();
    for (int iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); iLevel++)
    {
        for (auto patch : *hierarchy.getPatchLevel(iLevel))
        {
            auto guardedGrid = modelView.guardedGrid(*patch);
            checkIons(hybridModel->state.ions, getPatchPath(iLevel, *patch));
        }
    }
}

struct ElectromagHi5Diagnostic : public Hi5Diagnostic{
    ElectromagHi5Diagnostic() : Hi5Diagnostic{"electromag"}{}
};
TEST_F(ElectromagHi5Diagnostic, electromagWrite)
{
    dMan.addDiagDict(this->electromag("/EM_B")).addDiagDict(this->electromag("/EM_E")).dump();
}

struct ReadOnlyElectromagHi5Diagnostic : public Hi5Diagnostic{
    ReadOnlyElectromagHi5Diagnostic() : Hi5Diagnostic{"electromag", HighFive::File::ReadOnly}{}
};
TEST_F(ReadOnlyElectromagHi5Diagnostic, electromagRead)
{
    auto& hierarchy = basicHierarchy->getHierarchy();
    for (int iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); iLevel++)
    {
        for (auto patch : *hierarchy.getPatchLevel(iLevel))
        {
            auto guardedGrid = modelView.guardedGrid(*patch);
            std::string patchPath(getPatchPath(iLevel, *patch) + "/");
            auto& B = hybridModel->state.electromag.B;
            checkVecField(B, patchPath + B.name());
            auto& E = hybridModel->state.electromag.E;
            checkVecField(E, patchPath + E.name());
        }
    }
}

struct ParticlesHi5Diagnostic : public Hi5Diagnostic{
    ParticlesHi5Diagnostic() : Hi5Diagnostic{"particles"}{}
};
TEST_F(ParticlesHi5Diagnostic, particlesWrite)
{
    dMan.addDiagDict(this->particles("/ions/pop/ions_alpha/domain"))
        .addDiagDict(this->particles("/ions/pop/ions_alpha/levelGhost"))
        .addDiagDict(this->particles("/ions/pop/ions_alpha/patchGhost"))
        .addDiagDict(this->particles("/ions/pop/ions_protons/domain"))
        .addDiagDict(this->particles("/ions/pop/ions_protons/levelGhost"))
        .addDiagDict(this->particles("/ions/pop/ions_protons/patchGhost"))
        .dump();
}

struct ReadOnlyParticlesHi5Diagnostic : public Hi5Diagnostic{
    ReadOnlyParticlesHi5Diagnostic() : Hi5Diagnostic{"particles", HighFive::File::ReadOnly}{}
};
TEST_F(ReadOnlyParticlesHi5Diagnostic, particlesRead)
{
    auto checkParticles = [&](auto& particles, auto path) {
        if (!particles.size())
            return;
        std::vector<float> weightV, chargeV, vV;
        writer.file().getDataSet(path + "weight").read(weightV);
        writer.file().getDataSet(path + "charge").read(chargeV);
        writer.file().getDataSet(path + "v").read(vV);
        std::vector<int> iCellV;
        writer.file().getDataSet(path + "iCell").read(iCellV);
        std::vector<float> deltaV;
        writer.file().getDataSet(path + "delta").read(deltaV);

        ParticlePacker packer{particles};

        auto first       = packer.empty();
        size_t iCellSize = std::get<2>(first).size();
        size_t deltaSize = std::get<3>(first).size();
        size_t vSize     = std::get<4>(first).size();
        size_t part_idx  = 0;
        while (packer.hasNext())
        {
            auto next = packer.next();

            for (size_t i = 0; i < iCellSize; i++)
                EXPECT_EQ(iCellV[(part_idx * iCellSize) + i], std::get<2>(next)[i]);

            for (size_t i = 0; i < deltaSize; i++)
                EXPECT_FLOAT_EQ(deltaV[(part_idx * deltaSize) + i], std::get<3>(next)[i]);

            for (size_t i = 0; i < vSize; i++)
                EXPECT_FLOAT_EQ(vV[(part_idx * vSize) + i], std::get<4>(next)[i]);

            part_idx++;
        }
    };
    auto& hierarchy = basicHierarchy->getHierarchy();
    for (int iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); iLevel++)
    {
        for (auto patch : *hierarchy.getPatchLevel(iLevel))
        {
            auto guardedGrid = modelView.guardedGrid(*patch);
            std::string patchPath(getPatchPath(iLevel, *patch) + "/");
            for (auto& pop : hybridModel->state.ions)
            {
                std::string particlePath(patchPath + "/ions/pop/" + pop.name() + "/");
                checkParticles(pop.domainParticles(), particlePath + "domain/");
                checkParticles(pop.levelGhostParticles(), particlePath + "levelGhost/");
                checkParticles(pop.patchGhostParticles(), particlePath + "patchGhost/");
            }
        }
    }
}

TEST_F(ReadOnlyParticlesHi5Diagnostic, attributesRead)
{

    auto& hierarchy = basicHierarchy->getHierarchy();
    for (int iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); iLevel++)
    {
        for (auto patch : *hierarchy.getPatchLevel(iLevel))
        {
            auto guardedGrid = modelView.guardedGrid(*patch);
            std::string patchPath(getPatchPath(iLevel, *patch) + "/"), origin;
            writer.file().getGroup(patchPath).getAttribute("origin").read(origin);

            EXPECT_EQ(guardedGrid.grid_.origin().str(), origin);
        }
    }
}

int main(int argc, char** argv)
{
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
