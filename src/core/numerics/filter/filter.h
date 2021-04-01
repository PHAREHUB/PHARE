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

        // this is hard coded but probably will never change
        auto constexpr alpha = 0.5;

        // this is doing a memory allocation
        // which is not optimal since that function will be called
        // several time per patch. Another way would be that
        // the filter is owned by the Solver and has an interface
        // used by the solver to declare a resource to the resource manager
        // for that temporary object. This way the allocation would be done
        // for all patches only once per regrid. This is something that should
        // be done eventually.
        std::vector<typename Field::type> tmp{density.data(), density.data()+density.size()};

        if constexpr (Field::dimension == 1)
        {
            for (auto ix=start; ix <= end; ++ix)
            {
                // also the fact that 'tmp' is a vector (and not a field)
                // forces us to use [] and not (), which is ok in 1D but much less
                // extensible to nD than if we had one tmp Field registered.
                density(ix) = alpha*tmp[ix] + 0.5*(1-alpha)*(tmp[ix-1]+ tmp[ix+1]);
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
