
// input/input_1d_ratio_2.txt is unused but a reference

#include <mpi.h>

#include "include.h"

using namespace PHARE::diagnostic;
using namespace PHARE::diagnostic::h5;


template<typename Simulator>
struct FluidHi5Diagnostic : public Hi5Diagnostic<Simulator>
{
    FluidHi5Diagnostic(Simulator& simulator, unsigned flags)
        : Hi5Diagnostic<Simulator>{simulator, "fluid", flags}
    {
    }
};

TYPED_TEST(SimulatorTest, fluid)
{
    using GridLayout = typename TypeParam::PHARETypes::GridLayout_t;

    TypeParam sim;
    { // Scoped to destruct after dump
        FluidHi5Diagnostic<TypeParam> hi5{sim, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.fluid("/ions/density"))
            .addDiagDict(hi5.fluid("/ions/bulkVelocity"))
            .addDiagDict(hi5.fluid("/ions/pop/ions_alpha/density"))
            .addDiagDict(hi5.fluid("/ions/pop/ions_alpha/flux"))
            .addDiagDict(hi5.fluid("/ions/pop/ions_protons/density"))
            .addDiagDict(hi5.fluid("/ions/pop/ions_protons/flux"))
            .dump();
    }

    FluidHi5Diagnostic<TypeParam> hi5{sim, HighFive::File::ReadOnly};
    auto& writer      = hi5.writer;
    auto& hybridModel = *sim.getHybridModel();

    auto checkIons = [&](auto& ions, auto patchPath) {
        std::string path(patchPath + "/ions/");

        for (auto& pop : hi5.modelView.getIons())
        {
            std::string popPath(path + "pop/" + pop.name() + "/");
            hi5.checkField(pop.density(), popPath + "density");
            hi5.checkVecField(pop.flux(), popPath + "flux");
        }

        hi5.checkField(ions.density(), path + "density");
        hi5.checkVecField(ions.velocity(), path + "bulkVelocity");
    };

    auto visit = [&]([[maybe_unused]] GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        auto patchPath = writer.getPatchPath("time", iLevel, patchID);
        checkIons(hybridModel.state.ions, patchPath);
    };

    sim.visitHierarchy(visit, hybridModel);
}


template<typename Simulator>
struct ElectromagHi5Diagnostic : public Hi5Diagnostic<Simulator>
{
    ElectromagHi5Diagnostic(Simulator& simulator, unsigned flags)
        : Hi5Diagnostic<Simulator>{simulator, "electromag", flags}
    {
    }
};

TYPED_TEST(SimulatorTest, electromag)
{
    using GridLayout = typename TypeParam::PHARETypes::GridLayout_t;

    TypeParam sim;
    { // scoped to destruct after dump
        ElectromagHi5Diagnostic<TypeParam> hi5{sim, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.electromag("/EM_B")).addDiagDict(hi5.electromag("/EM_E")).dump();
    }

    ElectromagHi5Diagnostic<TypeParam> hi5{sim, HighFive::File::ReadOnly};
    auto& writer      = hi5.writer;
    auto& hybridModel = *sim.getHybridModel();

    auto visit = [&]([[maybe_unused]] GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        auto patchPath = writer.getPatchPath("time", iLevel, patchID) + "/";
        auto& B        = hybridModel.state.electromag.B;
        auto& E        = hybridModel.state.electromag.E;
        hi5.checkVecField(B, patchPath + B.name());
        hi5.checkVecField(E, patchPath + E.name());
    };

    sim.visitHierarchy(visit, hybridModel);
}


template<typename Simulator>
struct ParticlesHi5Diagnostic : public Hi5Diagnostic<Simulator>
{
    ParticlesHi5Diagnostic(Simulator& simulator, unsigned flags)
        : Hi5Diagnostic<Simulator>{simulator, "particles", flags}
    {
    }
};

TYPED_TEST(SimulatorTest, particles)
{
    using GridLayout = typename TypeParam::PHARETypes::GridLayout_t;

    TypeParam sim;
    { // scoped to destruct after dump
        ParticlesHi5Diagnostic<TypeParam> hi5{sim, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.particles("/ions/pop/ions_alpha/domain"))
            .addDiagDict(hi5.particles("/ions/pop/ions_alpha/levelGhost"))
            .addDiagDict(hi5.particles("/ions/pop/ions_alpha/patchGhost"))
            .addDiagDict(hi5.particles("/ions/pop/ions_protons/domain"))
            .addDiagDict(hi5.particles("/ions/pop/ions_protons/levelGhost"))
            .addDiagDict(hi5.particles("/ions/pop/ions_protons/patchGhost"))
            .dump();
    }

    ParticlesHi5Diagnostic<TypeParam> hi5{sim, HighFive::File::ReadOnly};
    auto& writer      = hi5.writer;
    auto& hybridModel = *sim.getHybridModel();

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

    auto visit = [&]([[maybe_unused]] GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        auto patchPath = writer.getPatchPath("time", iLevel, patchID) + "/";
        for (auto& pop : hybridModel.state.ions)
        {
            std::string particlePath(patchPath + "/ions/pop/" + pop.name() + "/");
            checkParticles(pop.domainParticles(), particlePath + "domain/");
            checkParticles(pop.levelGhostParticles(), particlePath + "levelGhost/");
        }
    };

    sim.visitHierarchy(visit, hybridModel);
}

TYPED_TEST(SimulatorTest, attributesRead)
{
    using GridLayout = typename TypeParam::PHARETypes::GridLayout_t;

    TypeParam sim;
    ParticlesHi5Diagnostic<TypeParam> hi5{sim, HighFive::File::ReadOnly};
    auto& writer      = hi5.writer;
    auto& hybridModel = *sim.getHybridModel();

    auto visit = [&](GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        std::string patchPath = writer.getPatchPath("time", iLevel, patchID) + "/", origin;
        writer.file().getGroup(patchPath).getAttribute("origin").read(origin);
        EXPECT_EQ(gridLayout.origin().str(), origin);
    };

    sim.visitHierarchy(visit, hybridModel);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
