
#ifndef PHARE_PHARE_INCLUDE_H
#define PHARE_PHARE_INCLUDE_H

#include "simulator/simulator.h"
#include "core/utilities/algorithm.h"

#include <iostream>

class SamraiLifeCycle
{
public:
    SamraiLifeCycle(int argc, char** argv)
    {
        SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();
    }
    ~SamraiLifeCycle()
    {
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
    }
};

struct RuntimeDiagnosticInterface
{
    RuntimeDiagnosticInterface(PHARE::ISimulator& _simulator)
        : simulator{_simulator}
    {
        auto dict        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
        auto dim         = dict["simulation"]["dimension"].template to<int>();
        auto interpOrder = dict["simulation"]["interp_order"].template to<int>();
        if (!PHARE::core::makeAtRuntime<Maker>(dim, interpOrder, Maker{*this}))
            throw std::runtime_error("Runtime diagnostic deduction failed");
    }

    struct Maker
    {
        Maker(RuntimeDiagnosticInterface& _rdi)
            : rdi{_rdi}
        {
        }

        template<typename Dimension, typename InterpOrder>
        bool operator()(std::size_t userDim, std::size_t userInterpOrder, Dimension dimension,
                        InterpOrder interp_order)
        {
            auto dict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
            if (dict["simulation"].contains("diagnostics"))
            {
                if (userDim == dimension() and userInterpOrder == interp_order())
                {
                    std::size_t constexpr d  = dimension();
                    std::size_t constexpr io = interp_order();

                    using PHARE_Types         = PHARE::PHARE_Types<d, io>;
                    using DiagnosticModelView = typename PHARE_Types::DiagnosticModelView;
                    using DiagnosticWriter    = typename PHARE_Types::DiagnosticWriter;

                    auto* simulator = dynamic_cast<PHARE::Simulator<d, io>*>(&rdi.simulator);

                    rdi.modelView = std::make_unique<DiagnosticModelView>(
                        *simulator->getPrivateHierarchy(), *simulator->getHybridModel());

                    rdi.writer = DiagnosticWriter::from(
                        *static_cast<DiagnosticModelView*>(rdi.modelView.get()),
                        dict["simulation"]["diagnostics"]);

                    rdi.dMan = PHARE::diagnostic::DiagnosticsManager<DiagnosticWriter>::from(
                        *static_cast<DiagnosticWriter*>(rdi.writer.get()),
                        dict["simulation"]["diagnostics"]);

                    return 1;
                }
            }
            return 0;
        }

        RuntimeDiagnosticInterface& rdi;
    };

    void dump() { dMan->dump(); }

    PHARE::ISimulator& simulator;
    std::unique_ptr<PHARE::diagnostic::IDiagnosticModelView> modelView;
    std::unique_ptr<PHARE::diagnostic::IDiagnosticWriter> writer;
    std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager> dMan;
};


#endif /*PHARE_PHARE_INCLUDE_H*/
