#ifndef PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_HPP

#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/partitionner/partitionner.hpp"

#include <cstddef>

namespace PHARE::core
{
template<auto alloc_mode, typename ParticleArray, std::size_t impl = 0>
class ParticleArrayPartitioner;
} // namespace PHARE::core

#include "particles_partitioning_cpu.hpp"

#if PHARE_HAVE_GPU

#include "particles_partitioning_gpu.hpp"

#else // no impl (yet)

namespace PHARE::core
{

template<typename ParticleArray> // force compile error showing no impl
class ParticleArrayPartitioner<AllocatorMode::GPU_UNIFIED, ParticleArray, /*impl = */ 0>;

} // namespace PHARE::core

#endif

#endif /*PHARE_CORE_DATA_PARTICLES_PARTITIONING_PARTICLES_PARTITIONING_HPP*/
