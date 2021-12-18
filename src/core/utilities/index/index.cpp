
#include "index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/types.hpp"

namespace PHARE
{
namespace core
{
    MeshIndex<1> make_index(std::uint32_t i) { return MeshIndex<1>(i); }

    MeshIndex<2> make_index(std::uint32_t i, std::uint32_t j) { return MeshIndex<2>(i, j); }

    MeshIndex<3> make_index(std::uint32_t i, std::uint32_t j, std::uint32_t k)
    {
        return MeshIndex<3>(i, j, k);
    }

} // namespace core
} // namespace PHARE
