
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


#include "types/amr_types.h"
#include "detail/highfive.h"

struct Hi5Diagnostic : public ::testing::Test, public AfullHybridBasicHierarchy
{
    auto dict(std::string&& type)
    {
        PHARE::initializer::PHAREDict<1> dict;
        dict["diag"]["name"]            = type;
        dict["diag"]["type"]            = type;
        dict["diag"]["species"]         = std::string{"SPECIES_" + type};
        dict["diag"]["compute_every"]   = std::size_t{1};
        dict["diag"]["write_every"]     = std::size_t{1};
        dict["diag"]["start_iteration"] = std::size_t{0};
        dict["diag"]["end_iteration"]   = std::numeric_limits<std::size_t>::max();
        return dict;
    }

    /*use object address as uuid in case of parallelism*/
    std::string filename() const
    {
        std::stringstream ss;
        ss << "/tmp/hi_" << this << ".5";
        return ss.str();
    }
};

TEST_F(Hi5Diagnostic, hdf5Electromag)
{
    using DiagnosticModelView
        = PHARE::AMRDiagnosticModelView<PHARE::amr::SAMRAI_Types, HybridModelT>;

    using DiagnosticWriter = PHARE::HighFiveDiagnostic<DiagnosticModelView>;

    DiagnosticModelView modelView{basicHierarchy->getHierarchy(), *hybridModel};
    DiagnosticWriter writer{modelView, filename()};
    PHARE::DiagnosticsManager<DiagnosticWriter> dMan{writer};
    dMan.addDiagDict(this->dict("electromag")).dump();

    auto checkField = [&](auto& vecField, auto vecFieldID, auto& path) {
        for (auto const key : {"x", "y", "z"})
        {
            auto& field = vecField.getComponent(PHARE::core::Components::at(key));
            std::string fieldPath(path + "/" + field.name());
            std::vector<double> readData;
            writer.file().getDataSet(fieldPath).read(readData);
            EXPECT_EQ(readData.size(), field.size());
            for (size_t i = 0; i < field.size(); i++)
                EXPECT_DOUBLE_EQ(readData[i], field.data()[i]);
        }
    };

    size_t patch_idx = 0;
    for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
    {
        auto guardedGrid = modelView.guardedGrid(*patch);
        std::stringstream patchID;
        patchID << patch->getGlobalId();
        std::string patchPath("/t#/pl0/p" + patchID.str());
        checkField(hybridModel->state.electromag.B, "B", patchPath);
        checkField(hybridModel->state.electromag.E, "E", patchPath);
        patch_idx++;
    }
}

TEST_F(Hi5Diagnostic, hdf5Particles)
{
    using DiagnosticModelView
        = PHARE::AMRDiagnosticModelView<PHARE::amr::SAMRAI_Types, HybridModelT>;
    using DiagnosticWriter = PHARE::HighFiveDiagnostic<DiagnosticModelView>;

    DiagnosticModelView modelView{basicHierarchy->getHierarchy(), *hybridModel};
    DiagnosticWriter writer{modelView, filename()};
    PHARE::DiagnosticsManager<DiagnosticWriter> dMan{writer};
    dMan.addDiagDict(this->dict("particles")).dump();

    auto checkParticle = [&](auto& particles, auto path) {
        if (!particles.size())
            return;
        std::vector<double> weightV, chargeV, vV;
        writer.file().getDataSet(path + "weight").read(weightV);
        writer.file().getDataSet(path + "charge").read(chargeV);
        writer.file().getDataSet(path + "v").read(vV);
        std::vector<int> iCellV;
        writer.file().getDataSet(path + "iCell").read(iCellV);
        std::vector<float> deltaV;
        writer.file().getDataSet(path + "delta").read(deltaV);

        PHARE::core::HasNextIterable<PHARE::ParticlePacker, PHARE::core::ParticleArray<1>> packer{
            particles};

        auto first       = packer.first();
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
                EXPECT_DOUBLE_EQ(vV[(part_idx * vSize) + i], std::get<4>(next)[i]);

            part_idx++;
        }
    };

    size_t patch_idx = 0;
    for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
    {
        auto guardedGrid = modelView.guardedGrid(*patch);
        std::stringstream patchID;
        patchID << patch->getGlobalId();
        std::string patchPath("/t#/pl0/p" + patchID.str());
        size_t pop_idx = 0;
        for (auto& pop : hybridModel->state.ions)
        {
            std::string particlePath(patchPath + "/ions/pop/" + pop.name() + "/");
            checkParticle(pop.domainParticles(), particlePath + "domain/");
            checkParticle(pop.levelGhostParticles(), particlePath + "lvlGhost/");
            checkParticle(pop.patchGhostParticles(), particlePath + "patchGhost/");
            pop_idx++;
        }
        patch_idx++;
    }
}

int main(int argc, char** argv)
{
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
