#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H

#include "detail/highfive.h"

namespace PHARE
{
namespace diagnostic
{
    namespace h5
    {
        template<typename HighFiveDiagnostic>
        class ElectromagDiagnosticWriter : public Hi5DiagnosticWriter<HighFiveDiagnostic>
        {
        public:
            using Hi5DiagnosticWriter<HighFiveDiagnostic>::hi5_;
            using Attributes = typename Hi5DiagnosticWriter<HighFiveDiagnostic>::Attributes;
            ElectromagDiagnosticWriter(HighFiveDiagnostic& hi5)
                : Hi5DiagnosticWriter<HighFiveDiagnostic>(hi5)
            {
            }
            void write(DiagnosticDAO&) override;
            void compute(DiagnosticDAO&) override {}
            void getDataSetInfo(std::string const& patchID, Attributes& patchAttributes) override;
            void initDataSets(std::vector<std::string> const& patchIDs,
                              Attributes& patchAttributes) override;

        private:
            size_t level_                                 = 0;
            static constexpr size_t num_electromag_fields = 2; // B and E
        };

        template<typename HighFiveDiagnostic>
        void
        ElectromagDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(std::string const& patchID,
                                                                       Attributes& patchAttributes)
        {
            auto& hi5     = this->hi5_;
            size_t fields = 0;
            for (auto* vecField : hi5.modelView().getElectromagFields())
            {
                std::string fieldID                       = "electromag_" + fields;
                patchAttributes[patchID][fieldID]["name"] = vecField->name();

                for (auto& [id, type] : core::Components::componentMap)
                {
                    patchAttributes[patchID][fieldID][id] = vecField->getComponent(type).size();
                }
                fields++;
            }

            patchAttributes[patchID]["electromag"]["fields"] = fields;
        }

        template<typename HighFiveDiagnostic>
        void ElectromagDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
            std::vector<std::string> const& patchIDs, Attributes& patchAttributes)
        {
            auto& hi5 = this->hi5_;

            auto initNullPatch = [&]() {
                for (size_t i = 0; i < num_electromag_fields; i++)
                    for (auto& [id, type] : core::Components::componentMap)
                        hi5.template createDataSet<float>("", 0);
            };

            auto initPatch = [&](auto& patchID, auto& attributes) {
                std::string path{hi5.getPatchPath("time", level_, patchID)};
                for (size_t i = 0; i < num_electromag_fields; i++)
                {
                    std::string fieldID = "electromag_" + i;

                    for (auto& [id, type] : core::Components::componentMap)
                    {
                        hi5.template createDataSet<float>(
                            path + "/"
                                + patchAttributes[patchID][fieldID]["name"]
                                      .template to<std::string>()
                                + "/" + id,
                            patchAttributes[patchID][fieldID][id].template to<size_t>());
                    }
                }
            };

            size_t patches    = hi5.patchIDs().size();
            size_t maxPatches = hi5.getMaxOfPerMPI(patches);
            for (size_t i = 0; i < patches; i++)
                initPatch(patchIDs[i], patchAttributes[patchIDs[i]]);
            for (size_t i = patches; i < maxPatches; i++)
                initNullPatch();
        }

        template<typename HighFiveDiagnostic>
        void ElectromagDiagnosticWriter<HighFiveDiagnostic>::write([
            [maybe_unused]] DiagnosticDAO& diagnostic)
        {
            auto& hi5 = this->hi5_;

            for (auto* vecField : hi5.modelView().getElectromagFields())
            {
                hi5.writeVecFieldAsDataset(hi5.patchPath() + "/" + vecField->name(), *vecField);
            }
        }
    } // namespace h5
} // namespace diagnostic
} // namespace PHARE

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H */
