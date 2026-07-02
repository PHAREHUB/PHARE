#ifndef PHARE_DIAGNOSTIC_MANAGER_HPP_
#define PHARE_DIAGNOSTIC_MANAGER_HPP_

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"

#include "initializer/data_provider.hpp"

#include "diagnostic_props.hpp"

#include <map>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>

namespace PHARE::diagnostic
{
enum class Mode { LIGHT, FULL };



template<typename DiagManager>
void registerDiagnostics(DiagManager& dMan, initializer::PHAREDict const& diagsParams)
{
    auto const diagTypes = []() {
        using Model = DiagManager::Model_t;

        if constexpr (solver::is_hybrid_model_v<Model>)
            return std::vector<std::string>{"fluid", "electromag", "particle", "meta", "info"};

        else if constexpr (solver::is_mhd_model_v<Model>)
            return std::vector<std::string>{"mhd", "meta", "electromag"};

        else
            static_assert(core::dependent_false_v<Model>, "Unsupported model type");

        throw std::runtime_error("registerDiagnostics: Shouldn't happen!");
    }();

    for (auto& diagType : diagTypes)
    {
        // several diags of the same type can be registeroed
        // corresponding to several blocks in user input
        // fluid0, fluid1, electromag0, electromag1, etc.
        std::size_t diagBlockID = 0;
        while (diagsParams.contains(diagType)
               && diagsParams[diagType].contains(diagType + std::to_string(diagBlockID)))
        {
            std::string const diagName = diagType + std::to_string(diagBlockID);
            dMan.addDiagDict(diagsParams[diagType][diagName]);
            ++diagBlockID;
        }
    }
}




class IDiagnosticsManager
{
public:
    virtual bool dump(double timeStamp, double timeStep)         = 0;
    virtual void dump_level(std::size_t level, double timeStamp) = 0;
    inline virtual ~IDiagnosticsManager();
};
IDiagnosticsManager::~IDiagnosticsManager() {}



template<typename Writer>
class DiagnosticsManager : public IDiagnosticsManager
{
public:
    using Model_t = typename Writer::Model_t;

    bool dump(double timeStamp, double timeStep) override;


    void dump_level(std::size_t level, double timeStamp) override;


    DiagnosticsManager(std::unique_ptr<Writer>&& writer_ptr)
        : writer_{std::move(writer_ptr)}
    {
        if (!writer_)
            throw std::runtime_error("Error: DiagnosticsManager received null Writer");
    }


    template<typename Hierarchy, typename Model>
    NO_DISCARD static std::unique_ptr<DiagnosticsManager>
    make_unique(Hierarchy& hier, Model& model, initializer::PHAREDict const& dict)
    {
        auto dMan = std::make_unique<DiagnosticsManager>(Writer::make_unique(hier, model, dict));
        registerDiagnostics(*dMan, dict);
        return dMan;
    }


    DiagnosticsManager& addDiagDict(initializer::PHAREDict const& dict);


    DiagnosticsManager& addDiagDict(initializer::PHAREDict&& dict) { return addDiagDict(dict); }


    NO_DISCARD auto& diagnostics() const { return diagnostics_; }


    NO_DISCARD Writer& writer() { return *writer_.get(); }


    DiagnosticsManager(DiagnosticsManager const&)            = delete;
    DiagnosticsManager(DiagnosticsManager&&)                 = delete;
    DiagnosticsManager& operator=(DiagnosticsManager const&) = delete;
    DiagnosticsManager& operator=(DiagnosticsManager&&)      = delete;

private:
    // A scheduled time is "reached" if it lies before the next step boundary: scheduledTime is
    // within (or behind) the current step [timeStamp, timeStamp+timeStep). Compared as
    // (scheduledTime - timeStamp) < timeStep with a float cast to truncate trailing fp imprecision
    // (same robustness as the historical |nextTime-timeStamp| < timeStep): a time exactly one step
    // ahead (scheduledTime == timeStamp+timeStep) is NOT consumed now -- it belongs to the next
    // step. There is no abs(): a time already behind timeStamp yields a large negative difference,
    // so it stays "reached" and the catch-up loop keeps advancing instead of freezing.
    NO_DISCARD bool reached_(double scheduledTime, double timeStamp, double timeStep) const
    {
        return static_cast<float>(scheduledTime - timeStamp) < static_cast<float>(timeStep);
    }

    // Advance idx past every scheduled time the current step has reached; return true if any.
    // The while-loop (vs a single ++) keeps the cadence from freezing when one step overshoots
    // several scheduled times (dt > period, or adaptive dt growing past the period): a single ++
    // would let nextTime fall more than a step behind currentTime, after which
    // |nextTime-currentTime| never drops below dt again and all remaining dumps are silently lost.
    NO_DISCARD bool catchUp_(std::vector<double> const& times, std::size_t& idx, double timeStamp,
                             double timeStep) const
    {
        bool acted = false;
        while (idx < times.size() and reached_(times[idx], timeStamp, timeStep))
        {
            acted = true;
            ++idx;
        }
        return acted;
    }

    bool needsElapsedAction_(double const nextTime) const
    {
        return core::mpi::unix_timestamp_now() > nextTime;
    }



    NO_DISCARD bool needsWrite_(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto const& diag_key     = diag.type + diag.quantity;
        auto& nextWriteTimestamp = nextWrite_[diag_key];
        auto& nextWriteElapsed   = nextWriteElapsed_[diag_key];

        auto const writeTimestampNow
            = catchUp_(diag.writeTimestamps, nextWriteTimestamp, timeStamp, timeStep);

        auto const writeElapsedNow
            = nextWriteElapsed < diag.elapsedTimestamps.size()
              and needsElapsedAction_(diag.elapsedTimestamps[nextWriteElapsed]);

        if (writeElapsedNow)
            ++nextWriteElapsed;

        return writeTimestampNow || writeElapsedNow;
    }


    NO_DISCARD bool needsCompute_(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        return catchUp_(diag.computeTimestamps, nextCompute_[diag.type + diag.quantity], timeStamp,
                        timeStep);
    }


    std::vector<DiagnosticProperties> diagnostics_;
    std::unique_ptr<Writer> writer_;
    std::map<std::string, std::size_t> nextCompute_;
    std::map<std::string, std::size_t> nextWrite_;
    std::map<std::string, std::size_t> nextWriteElapsed_;

    std::size_t iteration_ = 0; ///< coarse-step counter for writeNiterPeriod cadence

    std::time_t const start_time_{core::mpi::unix_timestamp_now()};
};



template<typename Writer>
DiagnosticsManager<Writer>&
DiagnosticsManager<Writer>::addDiagDict(initializer::PHAREDict const& diagParams)
{
    auto& diagProps           = diagnostics_.emplace_back(DiagnosticProperties{});
    diagProps.type            = diagParams["type"].template to<std::string>();
    diagProps.quantity        = diagParams["quantity"].template to<std::string>();
    diagProps.writeTimestamps = diagParams["write_timestamps"].template to<std::vector<double>>();

    if (diagParams.contains("elapsed_timestamps"))
    {
        diagProps.elapsedTimestamps
            = diagParams["elapsed_timestamps"].template to<std::vector<double>>();
        for (auto& time : diagProps.elapsedTimestamps)
            time += start_time_; // expected for comparison later
    }

    diagProps["flush_every"] = diagParams["flush_every"].template to<std::size_t>();

    if (diagParams.contains("write_niter_period"))
        diagProps.writeNiterPeriod = diagParams["write_niter_period"].template to<std::size_t>();

    diagProps.computeTimestamps
        = diagParams["compute_timestamps"].template to<std::vector<double>>();

    diagProps.nAttributes = diagParams["n_attributes"].template to<std::size_t>();
    for (std::size_t i = 0; i < diagProps.nAttributes; ++i)
    {
        std::string const idx = std::to_string(i);
        std::string const key = diagParams["attribute_" + idx + "_key"];
        diagProps.forward_file_attribute(key, diagParams["attribute_" + idx + "_value"]);
    }

    return *this;
}


template<typename Writer>
void DiagnosticsManager<Writer>::dump_level(std::size_t level, double timeStamp)
{
    std::vector<DiagnosticProperties*> activeDiagnostics;

    for (auto& diag : diagnostics_)
        activeDiagnostics.emplace_back(&diag);

    writer_->dump_level(level, activeDiagnostics, timeStamp);
}


template<typename Writer>
bool DiagnosticsManager<Writer>::dump(double timeStamp, double timeStep)
{
    std::vector<DiagnosticProperties*> activeDiagnostics;
    for (auto& diag : diagnostics_)
    {
        // iteration-based cadence: fires (compute + write) every writeNiterPeriod coarse steps.
        // Used by write_niter_period (the only timestamp-free option valid under adaptive dt).
        bool const niterNow
            = diag.writeNiterPeriod > 0 and (iteration_ % diag.writeNiterPeriod == 0);

        // needsCompute_ advances its own timestamp index (catch-up), so call it unconditionally
        bool const computeNow = needsCompute_(diag, timeStamp, timeStep);
        if (computeNow or niterNow)
            writer_->getDiagnosticWriterForType(diag.type)->compute(diag);
        // call needsWrite_ unconditionally so its timestamp index still advances
        bool const writeNow = needsWrite_(diag, timeStamp, timeStep);
        if (writeNow or niterNow)
        {
            activeDiagnostics.emplace_back(&diag);
        }
    }

    if (activeDiagnostics.size() > 0)
    {
        PHARE_LOG_SCOPE(1, "DiagnosticsManager::dump");
        writer_->dump(activeDiagnostics, timeStamp);
    }

    ++iteration_;
    return activeDiagnostics.size() > 0;
}

} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_MANAGER_HPP_ */
