
#ifndef PHARE_PHARE_INCLUDE_H
#define PHARE_PHARE_INCLUDE_H

#include "simulator/simulator.h"
#include "core/utilities/algorithm.h"

#include <iostream>

namespace PHARE
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(const std::string& message, const std::string& filename, const int line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};

class SamraiLifeCycle
{
public:
    SamraiLifeCycle(int argc = 0, char** argv = nullptr)
    {
        SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();

        std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
            = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
        SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
    }
    ~SamraiLifeCycle()
    {
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
    }

    static void reset()
    {
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::startup();
    }
};

class StaticSamraiLifeCycle : public SamraiLifeCycle
{
public:
    inline static StaticSamraiLifeCycle& INSTANCE()
    {
        static StaticSamraiLifeCycle i;
        return i;
    }
};

} // namespace PHARE

struct RuntimeDiagnosticInterface
{
    RuntimeDiagnosticInterface(PHARE::ISimulator& _simulator, PHARE::amr::Hierarchy& _hierarchy)
        : hierarchy{_hierarchy}
        , simulator{_simulator}
    {
        auto dict          = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
        auto dim           = dict["simulation"]["dimension"].template to<int>();
        auto interpOrder   = dict["simulation"]["interp_order"].template to<int>();
        auto nbRefinedPart = dict["simulation"]["refined_particle_nbr"].template to<int>();

        if (!PHARE::core::makeAtRuntime<Maker>(dim, interpOrder, nbRefinedPart, Maker{*this}))
            throw std::runtime_error("Runtime diagnostic deduction failed");
    }

    struct Maker
    {
        Maker(RuntimeDiagnosticInterface& _rdi)
            : rdi{_rdi}
        {
        }



        template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
        bool operator()(std::size_t userDim, std::size_t userInterpOrder, size_t userNbRefinedPart,
                        Dimension dimension, InterpOrder interp_order, NbRefinedPart nbRefinedPart)
        {
            auto dict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
            if (dict["simulation"].contains("diagnostics"))
            {
                if (userDim == dimension() and userInterpOrder == interp_order()
                    and userNbRefinedPart == nbRefinedPart())
                {
                    std::size_t constexpr d  = dimension();
                    std::size_t constexpr io = interp_order();
                    std::size_t constexpr nb = nbRefinedPart();

                    using PHARE_Types      = PHARE::PHARE_Types<d, io, nb>;
                    using ModelView        = typename PHARE_Types::ModelView;
                    using DiagnosticWriter = typename PHARE_Types::DiagnosticWriter;

                    auto* simulator = dynamic_cast<PHARE::Simulator<d, io, nb>*>(&rdi.simulator);
                    auto& hierarchy = rdi.hierarchy;

                    rdi.modelView
                        = std::make_unique<ModelView>(hierarchy, *simulator->getHybridModel());

                    rdi.writer
                        = DiagnosticWriter::from(*static_cast<ModelView*>(rdi.modelView.get()),
                                                 dict["simulation"]["diagnostics"]);

                    rdi.dMan = PHARE::diagnostic::DiagnosticsManager<DiagnosticWriter>::from(
                        *static_cast<DiagnosticWriter*>(rdi.writer.get()),
                        dict["simulation"]["diagnostics"]);

                    return true;
                }
            }
            return false;
        }

        RuntimeDiagnosticInterface& rdi;
    };

    void dump(double timestamp, double timestep) { dMan->dump(timestamp, timestep); }

    PHARE::amr::Hierarchy& hierarchy;
    PHARE::ISimulator& simulator;
    std::unique_ptr<PHARE::diagnostic::IModelView> modelView;
    std::unique_ptr<PHARE::diagnostic::IWriter> writer;
    std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager> dMan;
};


#endif /*PHARE_PHARE_INCLUDE_H*/
