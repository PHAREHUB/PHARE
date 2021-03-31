#ifndef FILTER_H
#define FILTER_H

#include "core/data/grid/gridlayoutdefs.h"

namespace PHARE::core {


class Filter
{
public:
    template <typename Field, typename GridLayout>
    void operator()(Field& density, GridLayout const& layout)
    {
        auto start = layout.physicalStartIndex(density, Direction::X);
        auto end = layout.physicalEndIndex(density, Direction::X);
        auto alpha = 0.5;
        if constexpr (Field::dimension == 1)
        {
            for (auto ix=start; ix <= end; ++ix)
            {
                density(ix) = alpha*density(ix) + 0.5*(1-alpha)*(density(ix-1)+ density(ix+1));
            }
        }
        else
        {
            std::runtime_error("unsupported filter dim");
        }

    }
};


}


#endif // FILTER_H
