
// input/input_1d_ratio_2.txt is unused but a reference

#include "test_diagnostics.hpp"

template<typename Simulator>
void fluid_test(Simulator&& sim, std::string out_dir)
{
    using HybridModel = typename Simulator::HybridModel;
    using Hierarchy   = typename Simulator::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.hierarchy;

    { // Scoped to destruct after dump
        Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, out_dir, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.fluid("/ions/charge_density"))
            .addDiagDict(hi5.fluid("/ions/bulkVelocity"))
            .addDiagDict(hi5.fluid("/ions/momentum_tensor"))
            .addDiagDict(hi5.fluid("/ions/pop/alpha/density"))
            .addDiagDict(hi5.fluid("/ions/pop/alpha/flux"))
            .addDiagDict(hi5.fluid("/ions/pop/alpha/momentum_tensor"))
            .addDiagDict(hi5.fluid("/ions/pop/protons/density"))
            .addDiagDict(hi5.fluid("/ions/pop/protons/flux"))
            .addDiagDict(hi5.fluid("/ions/pop/protons/momentum_tensor"));
        hi5.dump();
    }

    Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, out_dir,
                                              HighFive::File::ReadOnly};
    validateFluidDump(sim, hi5);
}

template<typename Simulator>
void electromag_test(Simulator&& sim, std::string out_dir)
{
    using HybridModel = typename Simulator::HybridModel;
    using Hierarchy   = typename Simulator::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.hierarchy;
    { // scoped to destruct after dump
        Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, out_dir, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.electromag("/EM_B")).addDiagDict(hi5.electromag("/EM_E"));
        hi5.dump();
    }

    Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, out_dir,
                                              HighFive::File::ReadOnly};
    validateElectromagDump(sim, hi5);
}


template<typename Simulator>
void particles_test(Simulator&& sim, std::string out_dir)
{
    using HybridModel = typename Simulator::HybridModel;
    using Hierarchy   = typename Simulator::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.hierarchy;

    { // scoped to destruct after dump
        Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, out_dir, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.particles("/ions/pop/alpha/domain"))
            .addDiagDict(hi5.particles("/ions/pop/alpha/levelGhost"))
            .addDiagDict(hi5.particles("/ions/pop/protons/domain"))
            .addDiagDict(hi5.particles("/ions/pop/protons/levelGhost"));
        hi5.dump();
    }

    Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, out_dir,
                                              HighFive::File::ReadOnly};
    validateParticleDump(sim, hi5);
}


template<typename Simulator>
void allFromPython_test(Simulator&& sim, std::string out_dir)
{
    using HybridModel = typename Simulator::HybridModel;
    using Hierarchy   = typename Simulator::Hierarchy;

    sim.dump(/*timestamp= */ 0, /*timestep= */ 1);
    // flush h5files, killing the DiagnosticsManager killes the H5TypeWriter shared_ptrs internally,
    //  and forces closed any open h5files, flushing any remaining contents, which we use below in
    //  this test
    sim.reset_dman();

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.hierarchy;

    Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, out_dir,
                                              HighFive::File::ReadOnly};

    validateFluidDump(sim, hi5);
    validateElectromagDump(sim, hi5);
    validateParticleDump(sim, hi5);
    validateAttributes(sim, hi5);
}
