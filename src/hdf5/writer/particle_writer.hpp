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
        auto constexpr dim            = Particles::dimension;
        using Packer                  = core::ParticlePacker<dim>;
        auto static const packer_keys = Packer::keys();

        auto write_ = [&](auto& soa) {
            h5file.template write_data_set_flat<2>(path + packer_keys[0], soa.weight().data());
            h5file.template write_data_set_flat<2>(path + packer_keys[1], soa.charge().data());
            h5file.template write_data_set(path + packer_keys[2], soa.iCell().data());
            h5file.template write_data_set(path + packer_keys[3], soa.delta().data());
            h5file.template write_data_set(path + packer_keys[4], soa.v().data());
        };

        if constexpr (Particles::is_contiguous)
        {
            write_(particles);
        }
        else
        {
            constexpr bool SOA = true;
            Packer packer(particles);
            core::ParticleArray<dim, SOA> copy{particles.size()};
            packer.pack(copy);
            write_(copy);
        }
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
