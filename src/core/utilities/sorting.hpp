#ifndef PHARE_SORTING_HPP
#define PHARE_SORTING_HPP

#include "core/utilities/box/box.hpp"

#include <cstddef>
#include <vector>

namespace PHARE::core
{

template<typename Array, std::size_t dim>
class CountingSort
{
public:
    CountingSort(std::size_t size = 0)
        : tmp_(size)
    {
    }

    using value_type = typename Array::value_type;
    using SortBox    = core::Box<int, dim>;

    void set_size(std::size_t size)
    {
        if (size > tmp_.size())
        {
            tmp_.resize(size);
        }
    }

    template<typename CellGetter>
    void sort(Array& toSort, CellGetter getCell)
    {
        // pay the price of reallocation if temporary
        // array is too small
        if (toSort.size() > tmp_.size())
        {
            set_size(toSort.size());
        }
        if (box_.isEmpty())
        {
            return;
        }

        // compute histogram
        for (int ip = 0; ip < tmp_.size(); ++ip)
        {
            auto const& item = toSort[ip];
            auto const& cell = getCell(item);
            hist_[cell]++;
        }

        int sum = 0;
        for (int i = 0; i < hist_.size(); ++i)
        {
            // all particles in cell i will be in [sum, sum+hist_[i])
            int const tmp = hist_[i];
            hist_[i]      = sum;
            sum += tmp;
        }

        for (int ip = 0; ip < toSort.size(); ++ip)
        {
            auto const& cell    = getCell(toSort[ip]);
            tmp_[hist_[cell]++] = toSort[ip];
        }

        // now put items back in toSort
        // in the right order
        for (int ip = 0; ip < toSort.size(); ++ip)
        {
            toSort[ip] = tmp_[ip];
        }
    }

private:
    std::vector<value_type> tmp_;
    SortBox box_;
    std::vector<int> hist_;
};




} // namespace PHARE::core



#endif
