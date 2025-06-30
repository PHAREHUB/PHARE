#ifndef PHARE_HDF5_PARTICLE_WRITER_HPP
#define PHARE_HDF5_PARTICLE_WRITER_HPP


#include "core/def.hpp"
#include "core/data/particles/particle_packer.hpp"

#include "hdf5/detail/hdf5_utils.hpp"
#include "hdf5/detail/h5/h5_file.hpp"

#include <vector>

namespace PHARE::hdf5
{
class ParticleWriter
{
public:
    template<typename H5File, typename Particles>
    static void write(H5File& h5file, Particles const& particles, std::string const& path)
    {
        constexpr auto dim              = Particles::dimension;
        using Packer                    = core::ParticlePacker<dim>;
        constexpr auto particle_members = Packer::empty();
        static auto& keys               = Packer::keys();


        Packer{particles}.pack_ranges_into([&](auto const& arr, auto const from) {
            auto const soa_members = arr();

            core::for_N<Packer::n_keys>([&](auto ki) {
                auto const [key, member] = std::get<ki>(soa_members);
                auto const actual        = std::get<ki>(particle_members);

                h5file.file()
                    .getDataSet(path + keys[ki])
                    .select({from, 0ul}, size_for<dim>(actual, arr.size()))
                    .write_raw(member.data());
            });
        });
    }



    template<std::size_t dim, typename T, typename Size>
    NO_DISCARD auto static size_for(T const& type, Size const& n_particles)
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
