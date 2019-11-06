
#include "kul/log.hpp"

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
<<<<<<< HEAD
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
=======
    ~Hi5Diagnostic() { KLOG(INF); }

>>>>>>> pre mpi merge for hdf5 tests
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
<<<<<<< HEAD

    std::string filename() const { return "hi5_test.5"; }

    std::string getPatchOrigin(GridYeeT& grid, HybridModelT& model)
    {
        auto& bieldx = model.state.electromag.B.getComponent(PHARE::core::Component::X);
        return grid
            .fieldNodeCoordinates(bieldx, grid.origin(),
                                  grid.ghostStartIndex(bieldx, PHARE::core::Direction::X))
            .str();
    }
<<<<<<< HEAD
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
=======

    std::string getPatchPath(GridYeeT& grid, HybridModelT& model)
    {
        return writer.getPatchPath("time", 0 /*level*/, getPatchOrigin(grid, model));
    }

    std::string getPatchPath(std::string globalId)
    {
        return writer.getPatchPath("time", 0 /*level*/, globalId);
    }

    template<typename Field>
    void checkField(Field& field, std::string path)
    {
        std::vector<double> fieldV;
        writer.file().getDataSet(path).read(fieldV);
        EXPECT_EQ(fieldV.size(), field.size());

        KLOG(INF) << path << " size" << fieldV.size();
        KLOG(INF) << path << " size" << field.size();
        KLOG(INF) << path << " hdf5[0] " << fieldV[0];
        KLOG(INF) << path << " orig[0] " << field.data()[0];

        for (size_t i = 0; i < fieldV.size(); i++)
            EXPECT_DOUBLE_EQ(fieldV[i], field.data()[i]);
    }
>>>>>>> pre mpi merge for hdf5 tests

    template<typename VecField>
    void checkVecField(VecField& vecField, std::string fieldPath)
    {
<<<<<<< HEAD
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
=======
        KLOG(INF) << fieldPath;
        for (auto& [id, type] : Components::componentMap)
        {
            checkField(vecField.getComponent(type), fieldPath + "/" + id);
        }
>>>>>>> pre mpi merge for hdf5 tests
    }

<<<<<<< HEAD
TEST_F(Hi5Diagnostic, hdf5Particles)
{
<<<<<<< HEAD
<<<<<<< HEAD
    using DiagnosticModelView
        = PHARE::AMRDiagnosticModelView<PHARE::amr::SAMRAI_Types, HybridModelT>;
    using DiagnosticWriter = PHARE::HighFiveDiagnostic<DiagnosticModelView>;
=======
=======
>>>>>>> pre mpi merge for hdf5 tests
    using DiagnosticModelView = PHARE::SamraiDiagnosticModelView<HybridModelT>;
    using DiagnosticWriter    = PHARE::HighFiveDiagnostic<DiagnosticModelView>;
>>>>>>> implement visitor pattern to remove samrai visibility for highfive diagnostics

    DiagnosticModelView modelView{basicHierarchy->getHierarchy(), *hybridModel};
    DiagnosticWriter writer{modelView, filename()};
    PHARE::DiagnosticsManager<DiagnosticWriter> dMan{writer};
};


TEST_F(Hi5Diagnostic, hdf5Electromag)
{
    try
    {
        dMan.addDiagDict(this->dict("electromag")).dump();

        // for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
        // {
        //     auto guardedGrid = modelView.guardedGrid(*patch);

<<<<<<< HEAD
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
=======
        //     std::stringstream patchID;
        //     patchID << patch->getGlobalId();
>>>>>>> pre mpi merge for hdf5 tests

        //     std::string path(getPatchPath(patchID.str()) + "/");
        //     auto& B = hybridModel->state.electromag.B;
        //     checkVecField(B, path + B.name());
        //     auto& E = hybridModel->state.electromag.E;
        //     checkVecField(E, path + E.name());
        // }
        KLOG(INF);
    }
    catch (std::exception& err)
    {
<<<<<<< HEAD
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
=======
        // catch and print any HDF5 error
        std::cerr << err.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
>>>>>>> pre mpi merge for hdf5 tests
    }
}

// TEST_F(Hi5Diagnostic, hdf5Particles)
// {
//     dMan.addDiagDict(this->dict("particles")).dump();

//     auto checkParticles = [&](auto& particles, auto path) {
//         if (!particles.size())
//             return;
//         std::vector<double> weightV, chargeV, vV;
//         writer.file().getDataSet(path + "weight").read(weightV);
//         writer.file().getDataSet(path + "charge").read(chargeV);
//         writer.file().getDataSet(path + "v").read(vV);
//         std::vector<int> iCellV;
//         writer.file().getDataSet(path + "iCell").read(iCellV);
//         std::vector<float> deltaV;
//         writer.file().getDataSet(path + "delta").read(deltaV);

//         PHARE::ParticlePacker packer{particles};

//         auto first       = packer.first();
//         size_t iCellSize = std::get<2>(first).size();
//         size_t deltaSize = std::get<3>(first).size();
//         size_t vSize     = std::get<4>(first).size();
//         size_t part_idx  = 0;
//         while (packer.hasNext())
//         {
//             auto next = packer.next();

//             for (size_t i = 0; i < iCellSize; i++)
//                 EXPECT_EQ(iCellV[(part_idx * iCellSize) + i], std::get<2>(next)[i]);

//             for (size_t i = 0; i < deltaSize; i++)
//                 EXPECT_FLOAT_EQ(deltaV[(part_idx * deltaSize) + i],
//                 std::get<3>(next)[i]);

//             for (size_t i = 0; i < vSize; i++)
//                 EXPECT_DOUBLE_EQ(vV[(part_idx * vSize) + i], std::get<4>(next)[i]);

//             part_idx++;
//         }
//     };

//     for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
//     {
//         auto guardedGrid = modelView.guardedGrid(*patch);
//         std::string patchPath(getPatchPath(guardedGrid, *hybridModel));
//         for (auto& pop : hybridModel->state.ions)
//         {
//             std::string particlePath(patchPath + "/ions/pop/" + pop.name() + "/");
//             checkParticles(pop.domainParticles(), particlePath + "domain/");
//             checkParticles(pop.levelGhostParticles(), particlePath + "lvlGhost/");
//             checkParticles(pop.patchGhostParticles(), particlePath + "patchGhost/");
//         }
//     }
// }



// TEST_F(Hi5Diagnostic, hdf5Fluid)
// {
//     dMan.addDiagDict(this->dict("fluid")).dump();

//     auto checkIons = [&](auto& ions, auto patchPath) {
//         std::string path(patchPath + "/ions/");

//         for (auto& pop : modelView.getIons())
//         {
//             std::string popPath(path + "pop/" + pop.name() + "/");
//             checkField(pop.density(), popPath + "density");
//             checkVecField(pop.flux(), popPath + "flux");
//         }

//         checkField(ions.density(), path + "density");
//         checkVecField(ions.velocity(), path + "bulkVelocity");
//     };

//     for (auto patch : *basicHierarchy->getHierarchy().getPatchLevel(0))
//     {
//         auto guardedGrid = modelView.guardedGrid(*patch);
//         std::string patchPath(getPatchPath(guardedGrid, *hybridModel));
//         checkIons(hybridModel->state.ions, patchPath);
//     }
// }



int main(int argc, char** argv)
{
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


//
// simple example to write a dataset with Parallel HDF5 with MPI-IO
//
// The dataset is written from several MPI node in parallel
//
//
// int main(int argc, char** argv)
// {
//     int mpi_rank, mpi_size;

//     // initialize MPI
//     // MPI_Init(&argc, &argv);
//     PHARE_test::SamraiLifeCycle samsam(argc, argv);
//     MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
//     MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

//     using namespace HighFive;
//     try
//     {
//         // open a new file with the MPI IO driver for parallel Read/Write
//         File file(FILE_NAME, File::ReadWrite | File::Create | File::Truncate,
//                   MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

//         // we define the size of our dataset to
//         //  lines : total number of mpi_rank
//         //  columns : 2
//         std::vector<size_t> dims(2);
//         dims[0] = std::size_t(mpi_size);
//         dims[1] = 2;

//         // Create the dataset
//         DataSet dataset = file.createDataSet<double>(DATASET_NAME, DataSpace(dims));

//         // Each node want to write its own rank two time in
//         // its associated row
//         int data[1][2] = {{mpi_rank, mpi_rank}};

//         // write it to the associated mpi_rank
//         dataset.select({std::size_t(mpi_rank), 0}, {1, 2}).write(data);
//     }
//     catch (Exception& err)
//     {
//         // catch and print any HDF5 error
//         std::cerr << err.what() << std::endl;
//         MPI_Abort(MPI_COMM_WORLD, 1);
//     }

//     // MPI_Finalize();
//     return 0; // successfully terminated
// }
