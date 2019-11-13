#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H

#include "detail/highfive.h"
#include "detail/highfive_diag_writer.h"

namespace PHARE
{
namespace diagnostic
{
    namespace h5
    {
        template<typename HighFiveDiagnostic>
        class FluidDiagnosticWriter : public Hi5DiagnosticWriter<HighFiveDiagnostic>
        {
        public:
            using Hi5DiagnosticWriter<HighFiveDiagnostic>::hi5_;
            using Attributes = typename Hi5DiagnosticWriter<HighFiveDiagnostic>::Attributes;
            FluidDiagnosticWriter(HighFiveDiagnostic& hi5)
                : Hi5DiagnosticWriter<HighFiveDiagnostic>(hi5)
            {
            }
            void write(DiagnosticDAO&) override;
            void compute(DiagnosticDAO&) override {}
            void getDataSetInfo(std::string const& patchID, Attributes& patchAttributes) override;
            void initDataSets(std::vector<std::string> const& patchIDs,
                              Attributes& patchAttributes) override;

        private:
            size_t maxPops = 0, level_ = 0;
        };

        template<typename HighFiveDiagnostic>
        void FluidDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(std::string const& patchID,
                                                                       Attributes& patchAttributes)
        {
            auto& hi5   = this->hi5_;
            size_t pops = 0;
            auto& ions  = hi5.modelView().getIons();
            for (auto& pop : ions)
            {
                std::string popId                          = "fluid_pop_" + pops;
                patchAttributes[patchID][popId]["name"]    = pop.name();
                patchAttributes[patchID][popId]["density"] = pop.density().size();

                for (auto& [id, type] : core::Components::componentMap)
                    patchAttributes[patchID][popId]["flux"][id]
                        = pop.flux().getComponent(type).size();
                pops++;
            }
            maxPops                                    = pops;
            patchAttributes[patchID]["ion"]["density"] = ions.density().size();
            for (auto& [id, type] : core::Components::componentMap)
                patchAttributes[patchID]["ion"]["bulkVelocity"][id]
                    = ions.velocity().getComponent(type).size();
        }

        template<typename HighFiveDiagnostic>
        void FluidDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
            std::vector<std::string> const& patchIDs, Attributes& patchAttributes)
        {
            auto& hi5 = this->hi5_;

            auto initNullPatch = [&]() {
                for (size_t i = 0; i < maxPops; i++)
                    for (auto& [id, type] : core::Components::componentMap)
                        hi5.template createDataSet<float>("", 0);
            };

            auto initPatch = [&](auto& patchID, auto& attributes) {
                std::string path{hi5.getPatchPath("time", level_, patchID) + "/ions/"};
                for (size_t i = 0; i < maxPops; i++)
                {
                    std::string popId = "fluid_pop_" + i;
                    std::string popName
                        = patchAttributes[patchID][popId]["name"].template to<std::string>();
                    std::string popPath(path + "pop/" + popName + "/");
                    hi5.template createDataSet<float>(
                        popPath + "density",
                        patchAttributes[patchID][popId]["density"].template to<size_t>());
                    for (auto& [id, type] : core::Components::componentMap)
                        hi5.template createDataSet<float>(
                            popPath + "flux" + "/" + id,
                            patchAttributes[patchID][popId]["flux"][id].template to<size_t>());
                }
                hi5.template createDataSet<float>(
                    path + "density",
                    patchAttributes[patchID]["ion"]["density"].template to<size_t>());
                for (auto& [id, type] : core::Components::componentMap)
                    hi5.template createDataSet<float>(
                        path + "bulkVelocity" + "/" + id,
                        patchAttributes[patchID]["ion"]["bulkVelocity"][id].template to<size_t>());
            };

            size_t patches    = hi5.patchIDs().size();
            size_t maxPatches = hi5.getMaxOfPerMPI(patches);
            for (size_t i = 0; i < patches; i++)
                initPatch(patchIDs[i], patchAttributes[patchIDs[i]]);
            for (size_t i = patches; i < maxPatches; i++)
                initNullPatch();
        }

        template<typename HighFiveDiagnostic>
        void
        FluidDiagnosticWriter<HighFiveDiagnostic>::write([[maybe_unused]] DiagnosticDAO& diagnostic)
        {
            auto& hi5  = this->hi5_;
            auto& ions = hi5.modelView().getIons();
            std::string path{hi5.patchPath() + "/ions/"};

            for (auto& pop : hi5.modelView().getIons())
            {
                std::string popPath(path + "pop/" + pop.name() + "/");
                auto& density = pop.density();
                hi5.writeDataSet(popPath + "density", density.data());
                hi5.writeVecFieldAsDataset(popPath + "flux", pop.flux());
            }

            auto& density = ions.density();
            hi5.writeDataSet(path + "density", density.data());
            hi5.writeVecFieldAsDataset(path + "bulkVelocity", ions.velocity());
        }
    } // namespace h5
} // namespace diagnostic
} // namespace PHARE

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H */
