
#include "kul/log.hpp"

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

#include "detail/highfive.h"

struct Hi5Diagnostic : public ::testing::Test, public AfullHybridBasicHierarchy
{
    ~Hi5Diagnostic() {}

    auto dict(std::string&& type)
    {
        PHARE::initializer::PHAREDict<1> dict;
        dict["name"]            = type;
        dict["type"]            = type;
        dict["species"]         = std::string{"SPECIES_" + type};
        dict["compute_every"]   = std::size_t{1};
        dict["write_every"]     = std::size_t{1};
        dict["start_iteration"] = std::size_t{0};
        dict["end_iteration"]   = std::numeric_limits<std::size_t>::max();
        return dict;
    }

    std::string filename() const { return "hi5_test.5"; }

    std::string getPatchOrigin(GridYeeT& grid, HybridModelT& model)
    {
        auto& bieldx = model.state.electromag.B.getComponent(PHARE::core::Component::X);
        return grid
            .fieldNodeCoordinates(bieldx, grid.origin(),
                                  grid.ghostStartIndex(bieldx, PHARE::core::Direction::X))
            .str();
    }

    std::string getPatchPath(PHARE::amr::SAMRAI_Types::patch_t& patch)
    {
        std::stringstream globalId;
        globalId << patch.getGlobalId();
        return writer.getPatchPath("time", 0 /*level*/, globalId.str());
    }

    template<typename Field>
    void checkField(Field& field, std::string path)
    {
        std::vector<double> fieldV;
        writer.file().getDataSet(path).read(fieldV);
        EXPECT_EQ(fieldV.size(), field.size());

        for (size_t i = 0; i < fieldV.size(); i++)
            EXPECT_DOUBLE_EQ(fieldV[i], field.data()[i]);
    }

    template<typename VecField>
    void checkVecField(VecField& vecField, std::string fieldPath)
    {
        for (auto& [id, type] : Components::componentMap)
        {
            checkField(vecField.getComponent(type), fieldPath + "/" + id);
        }
    }

    using DiagnosticModelView
        = PHARE::AMRDiagnosticModelView<PHARE::amr::SAMRAI_Types, HybridModelT>;
    using DiagnosticWriter = PHARE::HighFiveDiagnostic<DiagnosticModelView>;

    DiagnosticModelView modelView{basicHierarchy->getHierarchy(), *hybridModel};
    DiagnosticWriter writer{modelView, filename()};
    PHARE::DiagnosticsManager<DiagnosticWriter> dMan{writer};
};


TEST_F(Hi5Diagnostic, electromag)
{
    dMan.addDiagDict(this->dict("electromag")).dump();

    for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
    {
        auto guardedGrid = modelView.guardedGrid(*patch);
        std::string patchPath(getPatchPath(*patch) + "/");
        auto& B = hybridModel->state.electromag.B;
        checkVecField(B, patchPath + B.name());
        auto& E = hybridModel->state.electromag.E;
        checkVecField(E, patchPath + E.name());
    }
}

TEST_F(Hi5Diagnostic, particles)
{
    dMan.addDiagDict(this->dict("particles")).dump();

    auto checkParticles = [&](auto& particles, auto path) {
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

        PHARE::ParticlePacker packer{particles};

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

    for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
    {
        auto guardedGrid = modelView.guardedGrid(*patch);
        std::string patchPath(getPatchPath(*patch) + "/");
        for (auto& pop : hybridModel->state.ions)
        {
            std::string particlePath(patchPath + "/ions/pop/" + pop.name() + "/");
            checkParticles(pop.domainParticles(), particlePath + "domain/");
            checkParticles(pop.levelGhostParticles(), particlePath + "lvlGhost/");
            checkParticles(pop.patchGhostParticles(), particlePath + "patchGhost/");
        }
    }
}

TEST_F(Hi5Diagnostic, fluid)
{
    dMan.addDiagDict(this->dict("fluid")).dump();

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

    for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
    {
        auto guardedGrid = modelView.guardedGrid(*patch);
        std::string patchPath(getPatchPath(*patch) + "/");
        checkIons(hybridModel->state.ions, patchPath);
    }
}

int main(int argc, char** argv)
{
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
