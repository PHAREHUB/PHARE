#ifndef PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H
#define PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H

//#include "kul/dbg.hpp"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

#include "diagnostic_manager.h"
#include "diagnostic_dao.h"
#include "highfive_diag_writer.h"
#include "types/electromag.h"
#include "types/particle.h"
#include "types/fluid.h"

namespace PHARE
{
namespace diagnostic
{
    namespace h5
    {
        struct HighFiveFile
        {
            HighFive::File file_;
            PHARE::diagnostic::Mode mode_ = PHARE::diagnostic::Mode::LIGHT;

            static auto createHighFiveFile(std::string const path)
            {
                return HighFive::File
                {
                    path,
                        HighFive::File::ReadWrite | HighFive::File::Create
                            | HighFive::File::Truncate
#if defined(H5_HAVE_PARALLEL)
                        ,
                        HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL)
#endif
                };
            }

            HighFiveFile(std::string const path,
                         PHARE::diagnostic::Mode mode = PHARE::diagnostic::Mode::LIGHT)
                : file_{createHighFiveFile(path)}
                , mode_{mode}
            {
            }
            ~HighFiveFile() {}

            HighFive::File& file() { return file_; }

            HighFiveFile(const HighFiveFile&)             = delete;
            HighFiveFile(const HighFiveFile&&)            = delete;
            HighFiveFile& operator&(const HighFiveFile&)  = delete;
            HighFiveFile& operator&(const HighFiveFile&&) = delete;
        };




        template<typename ModelView>
        class HighFiveDiagnostic
        {
        public:
            using This       = HighFiveDiagnostic<ModelView>;
            using GridLayout = typename ModelView::GridLayout;
            using Attributes = typename ModelView::Attributes;

            static constexpr auto dimension   = ModelView::dimension;
            static constexpr auto interpOrder = GridLayout::interp_order;

            HighFiveDiagnostic(ModelView& modelView, std::string const hifivePath)
                : hi5_{hifivePath}
                , modelView_{modelView}
            {
            }

            ~HighFiveDiagnostic() {}

            void dump(std::vector<DiagnosticDAO*> const&);

            HighFive::File& file() { return hi5_.file(); }

            template<typename String>
            auto getDiagnosticWriterForType(String& writer)
            {
                return writers.at(writer);
            }

            std::string getPatchPath([[maybe_unused]] std::string time, int iLevel,
                                     std::string globalCoords)
            {
                return "/t#/pl" + std::to_string(iLevel) + "/p" + globalCoords;
            }

            auto& modelView() const { return modelView_; }
            const auto& patchPath() const { return patchPath_; }

            template<typename Type>
            void createDataSet(std::string const& path, size_t size)
            {
                if constexpr (std::is_same_v<Type, double>) // force doubles for floats for storage
                    this->template createDatasetsPerMPI<float>(path, size);
                else
                    this->template createDatasetsPerMPI<Type>(path, size);
            }

            template<typename Array, typename String>
            void writeDataSet(String path, Array const* const array)
            {
                hi5_.file_.getDataSet(path).write(array);
            }

            template<typename String, typename Data>
            void writeAttribute(String path, std::string const& key, Data const& value)
            {
                hi5_.file_.getGroup(path)
                    .template createAttribute<Data>(key, HighFive::DataSpace::From(value))
                    .write(value);
            }

            template<typename Array>
            void writeDataSetPart(std::string path, size_t start, size_t size, Array const& array)
            {
                hi5_.file_.getDataSet(path).select({start}, {size}).write(array);
            }

            template<typename VecField>
            void writeVecFieldAsDataset(std::string path, VecField& vecField)
            {
                for (auto& [id, type] : core::Components::componentMap)
                    hi5_.file_.getDataSet(path + "/" + id)
                        .write(vecField.getComponent(type).data());
            }

            size_t getMaxOfPerMPI(size_t localSize);

            const auto& patchIDs() const { return patchIDs_; }

        private:
            HighFiveFile hi5_;
            ModelView& modelView_;
            std::string patchPath_; // is passed around as "virtual write()" has no parameters
            std::vector<std::string> patchIDs_;

            std::unordered_map<std::string, std::shared_ptr<Hi5DiagnosticWriter<This>>> writers{
                {"fluid", make_writer<FluidDiagnosticWriter<This>>()},
                {"electromag", make_writer<ElectromagDiagnosticWriter<This>>()},
                {"particles", make_writer<ParticlesDiagnosticWriter<This>>()}};

            template<typename Writer>
            std::shared_ptr<Hi5DiagnosticWriter<This>> make_writer()
            {
                return std::make_shared<Writer>(*this);
            }

            template<typename Type>
            void createDatasetsPerMPI(std::string path, size_t dataSetSize);


            template<typename Dict>
            void writeDict(Dict, std::string);

            HighFiveDiagnostic(const HighFiveDiagnostic&)             = delete;
            HighFiveDiagnostic(const HighFiveDiagnostic&&)            = delete;
            HighFiveDiagnostic& operator&(const HighFiveDiagnostic&)  = delete;
            HighFiveDiagnostic& operator&(const HighFiveDiagnostic&&) = delete;
        };

        /*TO DO
         * investigate level > 0 for MPI
         *  finalise HDF5 Path format
         */
        template<typename ModelView>
        void HighFiveDiagnostic<ModelView>::dump(std::vector<DiagnosticDAO*> const& diagnostics)
        {
            writeAttribute("/", "dim", dimension);
            writeAttribute("/", "interpOrder", interpOrder);

            /*TODO
              add time/iterations
            */
            {
                Attributes patchAttributes; // stores dataset info/size for synced MPI creation
                auto visitPatch = [&](GridLayout& gridLayout, std::string patchID, int iLevel) {
                    patchIDs_.emplace_back(patchID);

                    for (auto* diagnostic : diagnostics)
                    {
                        writers.at(diagnostic->type)->getDataSetInfo({patchID}, patchAttributes);
                    }
                };
                modelView().visitHierarchy(visitPatch);

                for (auto* diagnostic : diagnostics)
                {
                    writers.at(diagnostic->type)->initDataSets(patchIDs_, patchAttributes);
                }
            }
            auto visitPatch = [&](GridLayout& gridLayout, std::string patchID, int iLevel) {
                patchPath_ = getPatchPath("time", 0, patchID);

                for (auto* diagnostic : diagnostics)
                {
                    writers.at(diagnostic->type)->write(*diagnostic);
                }
                writeDict(modelView().getPatchAttributes(gridLayout), patchPath_);
            };

            modelView().visitHierarchy(visitPatch);
        }

        /*
         * turns a dict of std::map<std::string, T> to hdf5 attributes
         */
        template<typename ModelView>
        template<typename Dict>
        void HighFiveDiagnostic<ModelView>::writeDict(Dict dict, std::string path)
        {
            using dict_map_t = typename Dict::map_t;
            auto visitor     = [&](auto&& map) {
                using Map_t = std::decay_t<decltype(map)>;
                if constexpr (std::is_same_v<Map_t, dict_map_t>)
                    for (auto& pair : map)
                        std::visit(
                            [&](auto&& val) {
                                using Val = std::decay_t<decltype(val)>;
                                if constexpr (core::is_dict_leaf<Val, Dict>::value)
                                    writeAttribute(path, pair.first, val);
                                else
                                    throw std::runtime_error(
                                        std::string("Expecting Writable value got ")
                                        + typeid(Map_t).name());
                            },
                            pair.second.get()->data);
                else
                    /* static_assert fails without if ! constexpr all possible types
                     * regardless of what it actually is.
                     */
                    throw std::runtime_error(std::string("Expecting map<string, T> got ")
                                             + typeid(Map_t).name());
            };
            std::visit(visitor, dict.data);
        }

        template<typename ModelView>
        size_t HighFiveDiagnostic<ModelView>::getMaxOfPerMPI(size_t localSize)
        {
            int mpi_size;
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
            uint64_t perMPI[mpi_size];
            MPI_Allgather(&localSize, 1, MPI_UINT64_T, perMPI, 1, MPI_UINT64_T, MPI_COMM_WORLD);
            for (const auto size : perMPI)
                if (size > localSize)
                    localSize = size;
            return localSize;
        }


        /*
         * Communicate all dataset paths and sizes to all MPI process to allow each to create all
         * datasets independently.  This is a requirement of HDF5. all paths must be the same length
         */
        template<typename ModelView>
        template<typename Type>
        void HighFiveDiagnostic<ModelView>::createDatasetsPerMPI(std::string path,
                                                                 size_t dataSetSize)
        {
            int mpi_size;
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

            uint64_t pathSize = path.size();
            uint64_t pathSizes[mpi_size], datasetSizes[mpi_size];
            MPI_Allgather(&pathSize, 1, MPI_UINT64_T, pathSizes, 1, MPI_UINT64_T, MPI_COMM_WORLD);
            MPI_Allgather(&dataSetSize, 1, MPI_UINT64_T, datasetSizes, 1, MPI_UINT64_T,
                          MPI_COMM_WORLD);

            uint64_t maxPathSize = pathSize;
            for (auto const& size : pathSizes)
                if (size > maxPathSize)
                    maxPathSize = size;

            char chars[maxPathSize * mpi_size];
            MPI_Allgather(path.c_str(), path.size(), MPI_CHAR, &chars, maxPathSize, MPI_CHAR,
                          MPI_COMM_WORLD);

            for (uint64_t i = 0; i < mpi_size; i++)
            {
                if (!datasetSizes[i])
                    continue;
                std::string datasetPath{&chars[pathSize * i], pathSize};
                H5Easy::detail::createGroupsToDataSet(hi5_.file_, datasetPath);
                hi5_.file_.createDataSet<Type>(datasetPath, HighFive::DataSpace(datasetSizes[i]));
            }
        }
    } // namespace h5
} // namespace diagnostic
} /* namespace PHARE */

#endif /* PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H */
