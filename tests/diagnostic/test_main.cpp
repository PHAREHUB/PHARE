
#include <mpi.h>

<<<<<<< HEAD
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "include.h"
#include "samrai_lifecycle.h"
<<<<<<< HEAD

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
=======
#include "diagnostic/include.h"
#include "diagnostic/samrai_lifecycle.h"
=======
>>>>>>> implement visitor pattern to remove samrai visibility for highfive diagnostics

#include "func.h"
#include "defs.h"
using namespace PHARE_test::_1d;

#include "integrator.h"
#include "tag_strat.h"
#include "hierarchy.h"

#include "diagnostic/samrai_diagnostic.h"
#include "detail/highfive.h"


struct Hi5Diagnostic : public AfullHybridBasicHierarchy
{
    HighFive::File file{PHARE::hi5::Diagnostic::createHighFiveFile("/tmp/hi.5")};
    ~Hi5Diagnostic() {}

>>>>>>> Enable Hi5 MPI if HDF5_IS_PARALLEL
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
<<<<<<< HEAD

    /*use object address as uuid in case of parallelism*/
    std::string filename() const
    {
        std::stringstream ss;
        ss << "/tmp/hi_" << this << ".5";
        return ss.str();
    }
=======
>>>>>>> Enable Hi5 MPI if HDF5_IS_PARALLEL
};

TEST_F(Hi5Diagnostic, hdf5Electromag)
{
<<<<<<< HEAD
<<<<<<< HEAD
    using DiagnosticModelView
        = PHARE::AMRDiagnosticModelView<PHARE::amr::SAMRAI_Types, HybridModelT>;

    using DiagnosticWriter = PHARE::HighFiveDiagnostic<DiagnosticModelView>;

    DiagnosticModelView modelView{basicHierarchy->getHierarchy(), *hybridModel};
    DiagnosticWriter writer{modelView, filename()};
    PHARE::DiagnosticsManager<DiagnosticWriter> dMan{writer};
=======
    using Model         = HybridModelT;
    using Samraive      = PHARE::SamraiHighFiveDiagnostic<Model>;
    constexpr auto mode = PHARE::diagnostic::Mode::LIGHT;

    PHARE::hi5::Diagnostic hi5{this->file, mode};
    Samraive samhighfo{basicHierarchy->getHierarchy(), *hybridModel, hi5};
    PHARE::DiagnosticsManager<Samraive> dMan{samhighfo};
>>>>>>> Enable Hi5 MPI if HDF5_IS_PARALLEL
=======
    using DiagnosticModelView = PHARE::SamraiDiagnosticModelView<HybridModelT>;
    using DiagnosticWriter    = PHARE::HighFiveDiagnostic<DiagnosticModelView>;

    DiagnosticModelView modelView{basicHierarchy->getHierarchy(), *hybridModel};
    DiagnosticWriter writer{modelView, filename()};
    PHARE::DiagnosticsManager<DiagnosticWriter> dMan{writer};
>>>>>>> implement visitor pattern to remove samrai visibility for highfive diagnostics
    dMan.addDiagDict(this->dict("electromag")).dump();

    auto checkField = [&](auto& vecField, auto vecFieldID, auto& path) {
        for (auto const key : {"x", "y", "z"})
        {
            auto& field = vecField.getComponent(PHARE::core::Components::at(key));
<<<<<<< HEAD
            std::string fieldPath(path + "/" + field.name());
            std::vector<double> readData;
            writer.file().getDataSet(fieldPath).read(readData);
            EXPECT_EQ(readData.size(), field.size());
            for (size_t i = 0; i < field.size(); i++)
                EXPECT_DOUBLE_EQ(readData[i], field.data()[i]);
=======
            std::stringstream fieldPath;
            fieldPath << path.str() << "/" << vecFieldID << key;
            std::vector<double> readData(field.size());
            this->file.getDataSet(fieldPath.str()).read(readData);
            ASSERT_TRUE(readData.size() == field.size());
            for (size_t i = 0; i < field.size(); i++)
                ASSERT_TRUE(readData[i] == field.data()[i]);
>>>>>>> Enable Hi5 MPI if HDF5_IS_PARALLEL
        }
    };

    size_t patch_idx = 0;
    for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
    {
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> implement visitor pattern to remove samrai visibility for highfive diagnostics
        auto guardedGrid = modelView.guardedGrid(*patch);
        std::stringstream patchID;
        patchID << patch->getGlobalId();
        std::string patchPath("/t#/pl0/p" + patchID.str());
<<<<<<< HEAD
=======
        auto guardedGrid = samhighfo.modelView().guardedGrid(*patch);
        std::stringstream patchPath;
        patchPath << "/t#/pl0/p" << patch_idx;
>>>>>>> Enable Hi5 MPI if HDF5_IS_PARALLEL
=======
>>>>>>> implement visitor pattern to remove samrai visibility for highfive diagnostics
        checkField(hybridModel->state.electromag.B, "B", patchPath);
        checkField(hybridModel->state.electromag.E, "E", patchPath);
        patch_idx++;
    }
}

TEST_F(Hi5Diagnostic, hdf5Particles)
{
<<<<<<< HEAD
<<<<<<< HEAD
    using DiagnosticModelView
        = PHARE::AMRDiagnosticModelView<PHARE::amr::SAMRAI_Types, HybridModelT>;
    using DiagnosticWriter = PHARE::HighFiveDiagnostic<DiagnosticModelView>;
=======
    using DiagnosticModelView = PHARE::SamraiDiagnosticModelView<HybridModelT>;
    using DiagnosticWriter    = PHARE::HighFiveDiagnostic<DiagnosticModelView>;
>>>>>>> implement visitor pattern to remove samrai visibility for highfive diagnostics

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
=======
    using Model         = HybridModelT;
    using Samraive      = PHARE::SamraiHighFiveDiagnostic<Model>;
    constexpr auto mode = PHARE::diagnostic::Mode::LIGHT;

    PHARE::hi5::Diagnostic hi5{this->file, mode};
    Samraive samhighfo{basicHierarchy->getHierarchy(), *hybridModel, hi5};
    PHARE::DiagnosticsManager<Samraive> dMan{samhighfo};
    dMan.addDiagDict(this->dict("particles")).dump();

    auto checkAttribute = [&](auto key, auto& path, auto const& original) {
        std::decay_t<decltype(original)> copy;
        this->file.getGroup(path).getAttribute(key).read(copy);
        ASSERT_TRUE(original == copy);
    };
    auto checkDataset = [&](auto key, auto& path, auto const& original) {
        std::decay_t<decltype(original)> copy;
        this->file.getDataSet(path + key).read(copy);
        ASSERT_TRUE(original == copy);
    };

    auto checkParticle = [&](auto& particle, auto path) {
        checkAttribute("weight", path, particle.weight);
        checkAttribute("charge", path, particle.charge);

        checkAttribute("Ex", path, particle.Ex);
        checkAttribute("Ey", path, particle.Ey);
        checkAttribute("Ez", path, particle.Ez);

        checkAttribute("Bx", path, particle.Bx);
        checkAttribute("By", path, particle.By);
        checkAttribute("Bz", path, particle.Bz);

        checkDataset("iCell", path, particle.iCell);
        checkDataset("delta", path, particle.delta);
        checkDataset("v", path, particle.v);
>>>>>>> Enable Hi5 MPI if HDF5_IS_PARALLEL
    };

    size_t patch_idx = 0;
    for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
    {
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> implement visitor pattern to remove samrai visibility for highfive diagnostics
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
=======
        auto guardedGrid = samhighfo.modelView().guardedGrid(*patch);
        std::stringstream patchPath;
        patchPath << "/t#/pl0/p" << patch_idx;
        size_t pop_idx = 0;
        for (auto& pop : hybridModel->state.ions)
        {
            std::stringstream particlePath;
            particlePath << patchPath.str() << "/pop/" << pop_idx << "/";
            size_t particle_idx = 0;
            for (auto& particle : pop.domainParticles())
            {
                std::stringstream domainPath;
                domainPath << particlePath.str() << "domain/" << particle_idx << "/";
                checkParticle(particle, domainPath.str());
                particle_idx++;
            }
>>>>>>> Enable Hi5 MPI if HDF5_IS_PARALLEL
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
