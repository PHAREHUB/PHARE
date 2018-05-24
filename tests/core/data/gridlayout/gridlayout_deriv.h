#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/field/field.h"
#include "data/grid/gridlayout.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "hybrid/hybrid_quantities.h"


using namespace PHARE;

std::vector<double> read(std::string filename)
{
    std::ifstream readFile(filename);
    assert(readFile.is_open());
    std::vector<double> x;

    std::copy(std::istream_iterator<double>(readFile), std::istream_iterator<double>(),
              std::back_inserter(x));
    return x;
}



/*
class a1DDerivative : public ::testing::Test
{
protected:
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> By;

public:
    a1DDerivative()
        : By{"By", PHARE::HybridQuantity::Scalar::By, 50}
    {
    }
};
*/



#endif
