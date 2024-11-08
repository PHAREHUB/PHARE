#ifndef PHARE_CLUSTERS_HPP
#define PHARE_CLUSTERS_HPP

#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"

#include <iostream>
#include <array>
#include <tuple>
#include <string>

namespace PHARE::core
{
template<typename Cluster>
class ClusterSet
{
public:
    static auto constexpr dimension = Cluster::dimension;

    ClusterSet(Box<int, dimension> const& box, std::array<int, dimension> const& cluster_size)
        : box_{box}
        , cluster_size_{cluster_size}
        , shape_{[&]() {
            std::array<int, dimension> s;
            auto bs = box.shape();
            for (auto i = 0u; i < dimension; ++i)
            {
                auto const div = (bs[i] + cluster_size_[i] - 1) / cluster_size_[i];
                s[i]           = div;
            }
            return s;
        }()}
        , clusters_(product(shape_))
    {
        for (auto idim = 0u; idim < dimension; ++idim)
        {
            if (box_.shape()[idim] < cluster_size_[idim])
            {
                throw std::runtime_error("cluster size larger than box size in dimension "
                                         + std::to_string(idim));
            }
        }



        auto const size_me = [&](auto dim, auto idx) {
            if (idx == shape_[dim] - 1)
            {
                auto const remain = box_.shape()[dim] % cluster_size_[dim];
                return (remain == 0) ? cluster_size_[dim] : remain;
            }
            else
                return cluster_size_[dim];
        };

        for (auto ix = 0u; ix < shape_[0]; ++ix)
        {
            if constexpr (dimension == 1)
            {
                // -1 because upper is included
                clusters_[ix].lower[0] = box.lower[0] + ix * cluster_size_[0];
                clusters_[ix].upper[0] = clusters_[ix].lower[0] + size_me(0, ix) - 1;
            }
            else
            {
                for (auto iy = 0u; iy < shape_[1]; ++iy)
                {
                    if constexpr (dimension == 2)
                    {
                        auto const i          = ix * shape_[1] + iy;
                        clusters_[i].lower[0] = box.lower[0] + ix * cluster_size_[0];
                        clusters_[i].upper[0] = clusters_[i].lower[0] + size_me(0, ix) - 1;
                        clusters_[i].lower[1] = box.lower[1] + iy * cluster_size_[1];
                        clusters_[i].upper[1] = clusters_[i].lower[1] + size_me(1, iy) - 1;
                    }
                    else
                    {
                        for (auto iz = 0u; iz < shape_[2]; ++iz)
                        {
                            auto const i = ix * shape_[1] * shape_[2] + shape_[2] * iy + iz;
                            clusters_[i].lower[0] = box.lower[0] + ix * cluster_size_[0];
                            clusters_[i].upper[0] = clusters_[i].lower[0] + size_me(0, ix) - 1;
                            clusters_[i].lower[1] = box.lower[1] + iy * cluster_size_[1];
                            clusters_[i].upper[1] = clusters_[i].lower[1] + size_me(1, iy) - 1;
                            clusters_[i].lower[2] = box.lower[2] + iz * cluster_size_[2];
                            clusters_[i].upper[2] = clusters_[i].lower[2] + size_me(2, iz) - 1;
                        }
                    }
                }
            }
        }
    }


    auto export_overlaped_with(Box<int, dimension> const& box) const
    {
        std::vector<std::pair<bool, Cluster>> overlaped;
        for (auto const& cluster : clusters_)
        {
            auto overlap = box * Box<int, dimension>{cluster.lower, cluster.upper};
            if (overlap)
            {
                auto complete_overlap = (*overlap).size() == cluster.size();
                overlaped.emplace_back(complete_overlap, cluster);
            }
        }
        return overlaped;
    }

    auto shape() const { return shape_; }
    auto size() const { return clusters_.size(); }

    auto begin() { return clusters_.begin(); }
    auto begin() const { return clusters_.begin(); }

    auto end() { return clusters_.end(); }
    auto end() const { return clusters_.end(); }

    auto& operator[](std::size_t i) { return clusters_[i]; }
    auto const& operator[](std::size_t i) const { return clusters_[i]; }

private:
    Box<int, dimension> box_;
    std::array<int, dimension> cluster_size_;
    std::array<int, dimension> shape_;
    std::vector<Cluster> clusters_;
};
} // namespace PHARE::core


#endif
