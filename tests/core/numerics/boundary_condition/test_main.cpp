#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <vector>

#include "core/data/particles/particle_array.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/utilities/box/box.hpp"



using namespace PHARE::core;

class ABoundaryConditionWhereAllParticlesLeave : public ::testing::Test
{
public:
    ABoundaryConditionWhereAllParticlesLeave()
        : boundaryBoxes{bbox}
        , leavingParticles_(10)
    {
        bc.setBoundaryBoxes(boundaryBoxes);
        for (auto& part : leavingParticles_)
        {
            part.iCell[0] = 5;  // these particles are out...
            part.iCell[1] = -1; // and not through the boundarybox
        }
    }


protected:
    static constexpr int interpOrder = 1;
    BoundaryCondition<2, interpOrder> bc;
    Box<int, 2> bbox{Point<int, 2>{0, 10}, Point<int, 2>{20, 11}};
    std::vector<Box<int, 2>> boundaryBoxes;
    ParticleArray<Particle<2>> leavingParticles_;
};




TEST_F(ABoundaryConditionWhereAllParticlesLeave, removesOutgoingParticles)
{
    auto toDelete
        = bc.applyOutgoingParticleBC(std::begin(leavingParticles_), std::end(leavingParticles_));
    leavingParticles_.erase(toDelete, std::end(leavingParticles_));

    EXPECT_EQ(0, leavingParticles_.size());
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
