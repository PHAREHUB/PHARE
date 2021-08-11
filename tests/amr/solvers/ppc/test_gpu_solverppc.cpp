#include <algorithm>
namespace unsafe {
    template <typename T>
    constexpr auto begin(const T& item) {
        return reinterpret_cast<const char*>(&item);
    }
    template <typename T>
    constexpr auto end(const T& item) {
        return reinterpret_cast<const char*>(&item)+sizeof(T);
    }
    template <typename T>
    auto bitwise_equal(const T& lhs, const T& rhs) {
        return std::equal(begin(lhs), end(lhs), // will become constexpr with C++20
                          begin(rhs), end(rhs));
    }
}

#include <string>
#include <iostream>

#include "amr/solvers/solver_ppc.h"
#include "amr/solvers/gpu/solver_ppc.h"
#include "tests/simulator/per_test.h"
#include "tests/core/data/field/test_field.h"

static std::string const job_file = "job_1d";

template<typename PHARE_TYPES>
struct Setup
{
    using Solver        = typename PHARE_TYPES::SolverPPC_t;
    using GridLayout    = typename PHARE_TYPES::GridLayout_t;
    using Hierarchy_t   = typename PHARE_TYPES::hierarchy_t;
    using HybridModel_t = typename PHARE_TYPES::HybridModel_t;
    using HybridState_t = typename HybridModel_t::State_t;
    using Types         = PHARE_TYPES;

    static constexpr auto dim           = PHARE_TYPES::dimension;
    static constexpr auto interp        = PHARE_TYPES::interp_order;
    static auto constexpr nbRefineParts = PHARE_TYPES::nbRefinedPart;

    auto static fill_states(Setup& self)
    {
        std::vector<PHARE::solver::gpu_mkn::PatchState<GridLayout>> states;
        PHARE::amr::visitHierarchy<GridLayout>(
            self.hierarchy, *self.hybridModel.resourcesManager,
            [&](auto& gridLayout, std::string patchID, size_t) {
                states.emplace_back(gridLayout, self.state);
            },
            self.topLvl, self.topLvl + 1, self.hybridModel);
        return states;
    }

    Setup(std::string job_id)
        : sim{job_id}
    {
    }

    SimulatorTestParam<dim, interp, nbRefineParts> sim;
    Hierarchy_t& hierarchy{*sim.hierarchy};
    HybridModel_t& hybridModel{*sim.getHybridModel()};
    HybridState_t& state{hybridModel.state};
    int topLvl{hierarchy.getNumberOfLevels() - 1};

    std::vector<PHARE::solver::gpu_mkn::PatchState<GridLayout>> states{fill_states(*this)};
};


int test_offloader(){    
    using PHARE_Types = PHARE::PHARE_Types<1, 1, 2, 1>;
    using Solver      = typename Setup<PHARE_Types>::Solver;
    using GridLayout  = typename Setup<PHARE_Types>::GridLayout;
    auto constexpr dim = GridLayout::dimension;
   
    std::size_t err = 0; 
  {
    Setup<PHARE_Types> setup{job_file};
    auto& dict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    Solver solver{dict["simulation"]["algo"]};
    PHARE::solver::gpu_mkn::Offloader<Solver> offloader{solver, dict};
    PHARE::amr::visitHierarchy<GridLayout>(
        setup.hierarchy, *setup.hybridModel.resourcesManager,
        [&](auto& gridLayout, std::string patchID, size_t) {
            offloader.alloc(gridLayout, setup.state);
        }, 0, 1, setup.hybridModel);
    offloader._move0(0.0001);
    offloader.visit([&](auto& batch, auto& patch_state, auto& pop, auto pop_idx){
      auto span = batch.get(0);
      for(std::size_t i = 0 ; i < pop.size(); i++){
        auto& p0 = span[i]; // copy back
        auto& p1 = pop.data()[i];
        if(p0 == p1) ++err;
        for(std::size_t dimdex = 0; dimdex < dim; ++dimdex)
          if(p0.iCell[dimdex] < p1.iCell[dimdex] - 1 ||
             p0.iCell[dimdex] > p1.iCell[dimdex] + 1) throw std::runtime_error("FAIL");
      }
    });
    offloader._move1(.0001);
    offloader.visit([&](auto& batch, auto& patch_state, auto& pop, auto pop_idx){
      auto span = batch.get(0);
      for(std::size_t i = 0 ; i < pop.size(); i++){
        auto& p0 = span[i]; // copy back
        auto& p1 = pop.data()[i];
        if(p0 == p1) ++err;
        for(std::size_t dimdex = 0; dimdex < dim; ++dimdex)
          if(p0.iCell[dimdex] < p1.iCell[dimdex] - 1 ||
             p0.iCell[dimdex] > p1.iCell[dimdex] + 1) throw std::runtime_error("FAIL");
      }
    });
  }
    KLOG(ERR) << err;
    return (err > 0) ? 1 : 0;
}

int test_offloader_noop(){    
    using PHARE_Types = PHARE::PHARE_Types<1, 1, 2, 1>;
    using Solver      = typename Setup<PHARE_Types>::Solver;
    using GridLayout  = typename Setup<PHARE_Types>::GridLayout;
    auto constexpr dim = GridLayout::dimension;
   
    std::size_t err = 0; 
  {
    Setup<PHARE_Types> setup{job_file};
    auto& dict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    Solver solver{dict["simulation"]["algo"]};
    PHARE::solver::gpu_mkn::Offloader<Solver> offloader{solver, dict};
    PHARE::amr::visitHierarchy<GridLayout>(
        setup.hierarchy, *setup.hybridModel.resourcesManager,
        [&](auto& gridLayout, std::string patchID, size_t) {
            offloader.alloc(gridLayout, setup.state);
        }, 0, 1, setup.hybridModel);
    offloader.template _move0<true>(0);

    auto as_vec=[](auto& arr){ return std::vector<double>(arr.data(), arr.data() + arr.size()); };    

    offloader.visit([&](auto& batch, auto& patch_state, auto& pop, auto pop_idx){
      auto span = batch.get(0);
      for(std::size_t i = 0 ; i < pop.size(); i++)
        if(!(span[i] == pop.data()[i])) ++err;

      auto layout = (*std::get<2>(batch))()[0];
      if(!unsafe::bitwise_equal(patch_state.layout, layout)) err += 11;
      auto& etc = std::get<3>(batch);
      using View = PHARE::core::FieldView<dim, double, double*, PHARE::core::HybridQuantity::Scalar>;
       
      PHARE::core::test(layout, 
        View{patch_state.density[pop_idx], PHARE::core::HybridQuantity::Scalar::rho}, 
        etc.density(), PHARE::FieldDomainFilter{});

      PHARE::core::test(layout, 
        View{patch_state.electromag[0], PHARE::core::HybridQuantity::Scalar::Ex}, 
        etc.Ex(), PHARE::FieldDomainFilter{});
      PHARE::core::test(layout, 
        View{patch_state.electromag[1], PHARE::core::HybridQuantity::Scalar::Ey}, 
        etc.Ey(), PHARE::FieldDomainFilter{});
      PHARE::core::test(layout, 
        View{patch_state.electromag[2], PHARE::core::HybridQuantity::Scalar::Ez}, 
        etc.Ez(), PHARE::FieldDomainFilter{});
        

      PHARE::core::test(layout, 
        View{patch_state.electromag[3], PHARE::core::HybridQuantity::Scalar::Bx}, 
        etc.Bx(), PHARE::FieldDomainFilter{});
      PHARE::core::test(layout, 
        View{patch_state.electromag[4], PHARE::core::HybridQuantity::Scalar::By}, 
        etc.By(), PHARE::FieldDomainFilter{});
      PHARE::core::test(layout, 
        View{patch_state.electromag[5], PHARE::core::HybridQuantity::Scalar::Bz}, 
        etc.Bz(), PHARE::FieldDomainFilter{});
        

      PHARE::core::test(layout, 
        View{patch_state.flux[0 + (3 * pop_idx)], PHARE::core::HybridQuantity::Scalar::rho}, 
        etc.fluxX(), PHARE::FieldDomainFilter{});
      PHARE::core::test(layout, 
        View{patch_state.flux[1 + (3 * pop_idx)], PHARE::core::HybridQuantity::Scalar::rho}, 
        etc.fluxY(), PHARE::FieldDomainFilter{});
      PHARE::core::test(layout, 
        View{patch_state.flux[2 + (3 * pop_idx)], PHARE::core::HybridQuantity::Scalar::rho}, 
        etc.fluxZ(), PHARE::FieldDomainFilter{});
        
    });
  }
    KLOG(ERR) << err;
    return (err > 0) ? 1 : 0;
}

int test(){
    using cpu_PHARE_Types = PHARE::PHARE_Types<1, 1, 2>;
    using gpu_PHARE_Types = PHARE::PHARE_Types<1, 1, 2, 1>;
    Setup<cpu_PHARE_Types> setup0{job_file};
    Setup<gpu_PHARE_Types> setup1{job_file};
    
//     auto& dict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
//     typename cpu_PHARE_Types::SolverPPC_t cpu_solver{dict};
//     typename gpu_PHARE_Types::SolverPPC_t cpu_solver{dict};
    
    setup0.sim.advance(.00001);
    setup1.sim.advance(.00001);

    PHARE::SamraiLifeCycle::reset();
    return 0;    
}

int main(int argc, char* argv[]){    
    PHARE::SamraiLifeCycle samsam(argc, argv);    
    return test_offloader_noop() + test_offloader();
}
