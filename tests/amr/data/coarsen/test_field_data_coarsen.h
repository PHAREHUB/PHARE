#ifndef PHARE_TEST_FIELD_DATA_COARSEN_H
#define PHARE_TEST_FIELD_DATA_COARSEN_H

#include "data/coarsening/field_coarsen.h"
#include "data/field/field_data_coarsen.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE;

using GridYee1DO1 = GridLayoutImplYee<1, 1>;
using Field1D     = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;

TEST(FieldDataCoarsen, canBeConstructed)
{
    std::size_t constexpr dimension = 1;
    FieldDataCoarsen<GridYee1DO1, Field1D> fieldDatacoarsen("linearCoarsenField");
}


#endif
