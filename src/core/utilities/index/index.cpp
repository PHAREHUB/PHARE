
#include "index.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/utilities/types.h"

namespace PHARE
{
namespace core
{
    MeshIndex<1> make_index(uint32 i) { return MeshIndex<1>(i); }

    MeshIndex<2> make_index(uint32 i, uint32 j) { return MeshIndex<2>(i, j); }

    MeshIndex<3> make_index(uint32 i, uint32 j, uint32 k) { return MeshIndex<3>(i, j, k); }

} // namespace core
} // namespace PHARE
