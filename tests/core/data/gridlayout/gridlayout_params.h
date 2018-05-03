#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_PARAMS_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_PARAMS_H

#include "data/grid/gridlayout.h"
#include "gridlayout_base_params.h"
#include "gridlayout_utilities.h"
#include "utilities/point/point.h"


#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <string>
#include <type_traits>

namespace PHARE
{
/* ::std::ostream& operator<<(::std::ostream& os, GridLayoutTestParam<Layout::Yee, 1> const& param);
 */
/* void PrintTo(GridLayoutTestParam<Layout::Yee, 1> const& param, ::std::ostream& os); */


inline HybridQuantity::Scalar getQuantity(uint32 iQuantity)
{
    switch (iQuantity)
    {
        case 0: return HybridQuantity::Scalar::Bx;
        case 1: return HybridQuantity::Scalar::By;
        case 2: return HybridQuantity::Scalar::Bz;
        case 3: return HybridQuantity::Scalar::Ex;
        case 4: return HybridQuantity::Scalar::Ey;
        case 5: return HybridQuantity::Scalar::Ez;
        case 6: return HybridQuantity::Scalar::Jx;
        case 7: return HybridQuantity::Scalar::Jy;
        case 8: return HybridQuantity::Scalar::Jz;
        case 9: return HybridQuantity::Scalar::rho;
        case 10: return HybridQuantity::Scalar::Vx;
        case 11: return HybridQuantity::Scalar::Vy;
        case 12: return HybridQuantity::Scalar::Vz;
        case 13: return HybridQuantity::Scalar::P;
        default:
            throw std::runtime_error("There is no quantity indexed by "
                                     + std::to_string(iQuantity));
    }
}



template<Layout layout, std::size_t dim>
auto createParam(std::string const& layoutName, uint32 interpOrder,
                 std::array<double, dim> const& dxdydz, std::array<uint32, dim> const& nbCellXYZ,
                 Point<double, dim> const& origin)
{
    GridLayoutTestParam<layout, dim> param{};

    param.layout = std::make_shared<GridLayout<layout, dim>>(dxdydz, nbCellXYZ, layoutName, origin,
                                                             interpOrder);

    param.interpOrder = interpOrder;



    param.dxdydz    = dxdydz;
    param.nbCellXYZ = nbCellXYZ;
    param.origin    = origin;


    return param;
}


template<typename Array>
void writeToArray(std::ifstream& stream, Array& array)
{
    std::size_t dim = array.size();
    stream >> array[0];
    if (dim > 1)
    {
        stream >> array[1];
    }
    if (dim > 2)
    {
        stream >> array[2];
    }
}


} // namespace PHARE




#endif
