
#include <mpi.h>

#include "diagnostic/include.h"
#include "diagnostic/samrai_lifecycle.h"

#include "diagnostic/func.h"
#include "diagnostic/defs.h"
using namespace PHARE_test::_1d;

#include "diagnostic/detail/samrai_highfive.h"

#include "diagnostic/integrator.h"
#include "diagnostic/tag_strat.h"
#include "diagnostic/hierarchy.h"

struct Hi5Diagnostic : public AfullHybridBasicHierarchy
{
    HighFive::File file{PHARE::hi5::Diagnostic::createHighFiveFile(filename())};

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
    using Model         = HybridModelT;
    using Samraive      = PHARE::SamraiHighFiveDiagnostic<Model>;
    constexpr auto mode = PHARE::diagnostic::Mode::LIGHT;

    PHARE::hi5::Diagnostic hi5{this->file, mode};
    Samraive samhighfo{basicHierarchy->getHierarchy(), *hybridModel, hi5};
    PHARE::DiagnosticsManager<Samraive> dMan{samhighfo};
    dMan.addDiagDict(this->dict("electromag")).dump();

    auto checkField = [&](auto& vecField, auto vecFieldID, auto& path) {
        for (auto const key : {"x", "y", "z"})
        {
            auto& field = vecField.getComponent(PHARE::core::Components::at(key));
            std::stringstream fieldPath;
            fieldPath << path.str() << "/" << vecFieldID << key;
            std::vector<double> readData(field.size());
            this->file.getDataSet(fieldPath.str()).read(readData);
            ASSERT_TRUE(readData.size() == field.size());
            for (size_t i = 0; i < field.size(); i++)
                ASSERT_TRUE(readData[i] == field.data()[i]);
        }
    };

    size_t patch_idx = 0;
    for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
    {
        auto guardedGrid = samhighfo.modelView().guardedGrid(*patch);
        std::stringstream patchPath;
        patchPath << "/t#/pl0/p" << patch_idx;
        checkField(hybridModel->state.electromag.B, "B", patchPath);
        checkField(hybridModel->state.electromag.E, "E", patchPath);
        patch_idx++;
    }
}

TEST_F(Hi5Diagnostic, hdf5Particles)
{
    using Model         = HybridModelT;
    using Samraive      = PHARE::SamraiHighFiveDiagnostic<Model>;
    constexpr auto mode = PHARE::diagnostic::Mode::LIGHT;

    PHARE::hi5::Diagnostic hi5{this->file, mode};
    Samraive samhighfo{basicHierarchy->getHierarchy(), *hybridModel, hi5};
    PHARE::DiagnosticsManager<Samraive> dMan{samhighfo};
    dMan.addDiagDict(this->dict("particles")).dump();

    auto checkParticle = [&](auto& particles, auto path) {
        std::vector<double> weightV, chargeV, vV;
        this->file.getDataSet(path + "weight").read(weightV);
        this->file.getDataSet(path + "charge").read(chargeV);
        this->file.getDataSet(path + "v").read(vV);
        std::vector<int> iCellV;
        this->file.getDataSet(path + "iCell").read(iCellV);
        std::vector<float> deltaV;
        this->file.getDataSet(path + "delta").read(deltaV);

        PHARE::ParticularPacker<PHARE::ParticlePackerPart> packer(particles);

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
        auto guardedGrid = samhighfo.modelView().guardedGrid(*patch);
        std::stringstream patchPath;
        patchPath << "/t#/pl0/p" << patch_idx;
        size_t pop_idx = 0;
        for (auto& pop : hybridModel->state.ions)
        {
            std::stringstream particlePath;
            particlePath << patchPath.str() << "/ions/pop/" << pop_idx << "/";
            {
                std::stringstream domainPath;
                domainPath << particlePath.str() << "domain/";
                checkParticle(pop.domainParticles(), domainPath.str());
            }
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
