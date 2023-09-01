#ifndef PHARE_SORTING_HPP
#define PHARE_SORTING_HPP

#include "core/utilities/box/box.hpp"

#include <iostream>
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

    void setup(std::size_t nbr_elements, SortBox const& box)
    {
        tmp_.resize(nbr_elements);
        hist_.resize(box.size());
    }

    template<typename CellGetter>
    void sort(Array& toSort, CellGetter getCell)
    {
        assert(toSort.size() == tmp_.size());
        // compute histogram
        for (std::size_t ip = 0; ip < tmp_.size(); ++ip)
        {
            auto const& item = toSort[ip];
            auto const& cell = getCell(item);
            hist_[cell]++;
        }

        int sum = 0;
        for (std::size_t i = 0; i < hist_.size(); ++i)
        {
            // all particles in cell i will be in [sum, sum+hist_[i])
            int const tmp = hist_[i];
            hist_[i]      = sum;
            sum += tmp;
        }

        for (std::size_t ip = 0; ip < toSort.size(); ++ip)
        {
            auto const& cell    = getCell(toSort[ip]);
            tmp_[hist_[cell]++] = toSort[ip];
        }

        // now put items back in toSort
        // in the right order
        for (std::size_t ip = 0; ip < toSort.size(); ++ip)
        {
            toSort[ip] = tmp_[ip];
        }
    }

private:
    std::vector<value_type> tmp_;
    std::vector<int> hist_;
};




} // namespace PHARE::core



#endif
