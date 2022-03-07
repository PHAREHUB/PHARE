#ifndef PHARE_HDF5_PARTICLE_WRITER_HPP
#define PHARE_HDF5_PARTICLE_WRITER_HPP

#include <array>
#include <vector>

#include "hdf5/detail/hdf5_utils.hpp"
#include "hdf5/detail/h5/h5_file.hpp"

#include "core/data/particles/particle_packer.hpp"


namespace PHARE::hdf5
{
class ParticleWriter
{
public:
    template<typename H5File, typename Particles>
    static void write(H5File& h5file, Particles const& particles, std::string const& path)
    {
        auto constexpr dim = Particles::dimension;
        using Packer       = core::ParticlePacker<dim>;

        Packer packer(particles);
        core::ContiguousParticles<dim> copy{particles.size()};
        packer.pack(copy);

        std::size_t part_idx = 0;
        core::apply(copy.as_tuple(), [&](auto const& arg) {
            auto data_path = path + packer.keys()[part_idx++];
            h5file.template write_data_set_flat<2>(data_path, arg.data());
        });
    }




    template<std::size_t dim, typename T, typename Size>
    auto static size_for(T const& type, Size const& n_particles)
    {
        if (n_particles == 0)
            return std::vector<std::size_t>{0};
        if constexpr (is_array_dataset<T, dim>)
            return std::vector<std::size_t>{n_particles, type.size()};
        else /* not an array so value one of type T*/
            return std::vector<std::size_t>{n_particles, 1};
    }
};



} // namespace PHARE::hdf5


#endif /* PHARE_HDF5_PARTICLE_WRITER_HPP */
