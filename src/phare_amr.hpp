#ifndef PHARE_AMR_INCLUDE_HPP
#define PHARE_AMR_INCLUDE_HPP


#include "phare_mpi.hpp" // IWYU pragma: keep

#include "amr/wrappers/hierarchy.hpp"




namespace PHARE::amr
{
template<auto opts>
struct PHARE_Types
{
    auto static constexpr dimension = opts.dimension;

    using hierarchy_t = PHARE::amr::Hierarchy;
};

} // namespace PHARE::amr


#endif // PHARE_AMR_INCLUDE_HPP
