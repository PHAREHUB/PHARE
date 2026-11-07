#ifndef PHARE_AMR_DATA_PARTICLES_PARTICLES_RESTART_HPP
#define PHARE_AMR_DATA_PARTICLES_PARTICLES_RESTART_HPP

// particle arrays to and from SAMRAI restart databases without a copy: each Particle member is
//  a strided hyperslab of the AoS in memory, HDF5 gathers/scatters it to/from the per member
//  datasets SAMRAI's putDoubleArray/putIntegerArray would write ("<name>_weight", "_iCell",
//  "_delta", "_v", big endian, SAMRAI's "Type" attribute), so files are as they were. HDF5
//  stages through its type conversion buffer, 1MiB by default (H5Pset_buffer), not a copy of
//  the particles.
//
// SAMRAI restarts are HDFDatabase: RestartManager's only default factory, and PHARE has no
//  restarts without HighFive (NullOpRestartsManager), anything else throws.
//  Without HDF5 SAMRAI has no restart database factory


#include "core/data/particles/particle_packer.hpp"

// SAMRAI's HAVE_HDF5, not PHARE_HAS_HIGHFIVE: that is per target and libraries (phare_amr)
//  include this without it, ParticlesData would differ between translation units
#include <SAMRAI/SAMRAI_config.h>
#include <SAMRAI/tbox/Database.h>

#ifdef HAVE_HDF5
#include <SAMRAI/tbox/HDFDatabase.h>
#endif

#include <tuple>
#include <cstddef>
#include <string>
#include <stdexcept>
#include <type_traits>

namespace PHARE::amr
{
class ParticlesRestart
{
public:
    // nothing for no particles: SAMRAI errors on 0 size arrays, get() then finds no keys
    template<typename Particles>
    static void put(SAMRAI::tbox::Database& db, std::string const& name, Particles const& particles)
    {
        if (particles.size() == 0)
            return;

        auto const group = group_id_(db);
        for_each_member_<Particles>(
            [&](auto const& key, auto type, auto const offset, auto const N) {
                using T = decltype(type);
                put_member_<T, Particles>(group, name + "_" + key, particles, offset, N);
            });
    }

    // particles must be empty, stay empty if name has no keys, read in place and mapped
    template<typename Particles>
    static void get(SAMRAI::tbox::Database& db, std::string const& name, Particles& particles)
    {
        if (particles.size() != 0)
            throw std::runtime_error("ParticlesRestart::get: particles not empty");

        std::size_t n_keys = 0, n_found = 0;
        for_each_member_<Particles>([&](auto const& key, auto...) {
            ++n_keys;
            n_found += db.keyExists(name + "_" + key);
        });
        if (n_found == 0) // can't read what doesn't exist
            return;
        if (n_found != n_keys)
            throw std::runtime_error("ParticlesData::getFromRestart has been given an "
                                     "invalid input file, inconsistent state detected");

        auto const group = group_id_(db);
        particles.resize(db.getArraySize(
            name + "_" + core::ParticlePacker<Particles::dimension>::arbitrarySingleValueKey()));
        for_each_member_<Particles>(
            [&](auto const& key, auto type, auto const offset, auto const N) {
                using T = decltype(type);
                get_member_<T, Particles>(group, name + "_" + key, particles, offset, N);
            });
        particles.map_particles(); // resize bypasses the cell map
    }

private:
    // fn(key, value type{}, byte offset in Particle, values per particle) per member, as
    //  ParticlePacker::get/keys: the restart file format
    template<typename Particles>
    static void for_each_member_(auto&& fn)
    {
        using Particle_t = Particles::Particle_t;
        using Packer     = core::ParticlePacker<Particles::dimension>;
        static_assert(std::is_standard_layout_v<Particle_t>, "member offsets need standard layout");

        Particle_t const particle{};
        auto const base      = reinterpret_cast<std::byte const*>(&particle);
        std::size_t key      = 0;
        auto const member_fn = [&](auto const& member) {
            using Member = std::decay_t<decltype(member)>;
            auto const offset
                = static_cast<std::size_t>(reinterpret_cast<std::byte const*>(&member) - base);
            if constexpr (std::is_arithmetic_v<Member>)
                fn(Packer::keys()[key++], Member{}, offset, std::size_t{1});
            else
                fn(Packer::keys()[key++], typename Member::value_type{}, offset,
                   std::tuple_size_v<Member>);
        };
        std::apply([&](auto const&... members) { (member_fn(members), ...); },
                   Packer::get(particle));
    }

#ifdef HAVE_HDF5
    static hid_t group_id_(SAMRAI::tbox::Database& db)
    {
        if (auto* hdf = dynamic_cast<SAMRAI::tbox::HDFDatabase*>(&db))
            return hdf->getGroupId();
        throw std::runtime_error("ParticlesRestart: restart database is not HDFDatabase");
    }

    static void check_(auto const id, std::string const& what)
    {
        if (id < 0)
            throw std::runtime_error("ParticlesRestart: HDF5 " + what + " failed");
    }

    // N consecutive T of each particle, offset bytes into it
    template<typename T, typename Particles>
    static hid_t memspace_(std::size_t const n_particles, std::size_t const offset,
                           std::size_t const N)
    {
        using Particle_t = Particles::Particle_t;
        static_assert(sizeof(Particle_t) % sizeof(T) == 0);
        if (offset % sizeof(T) != 0)
            throw std::runtime_error("ParticlesRestart: misaligned particle member");

        hsize_t const size   = n_particles * (sizeof(Particle_t) / sizeof(T));
        hsize_t const start  = offset / sizeof(T);
        hsize_t const stride = sizeof(Particle_t) / sizeof(T);
        hsize_t const count  = n_particles;
        hsize_t const block  = N;

        hid_t const space = H5Screate_simple(1, &size, nullptr);
        check_(space, "H5Screate_simple");
        check_(H5Sselect_hyperslab(space, H5S_SELECT_SET, &start, &stride, &count, &block),
               "H5Sselect_hyperslab");
        return space;
    }

    // as SAMRAI's HDFDatabase: H5T_SAMRAI_DOUBLE / H5T_SAMRAI_INT
    template<typename T>
    static hid_t file_type_()
    {
        if constexpr (std::is_same_v<T, double>)
            return H5T_IEEE_F64BE;
        else
            return H5T_STD_I32BE;
    }

    template<typename T>
    static hid_t native_type_()
    {
        if constexpr (std::is_same_v<T, double>)
            return H5T_NATIVE_DOUBLE;
        else
            return H5T_NATIVE_INT;
    }

    // SAMRAI's HDFDatabase "Type" attribute, its private KEY_DOUBLE_ARRAY = 5 and
    //  KEY_INT_ARRAY = 7 (HDFDatabase.cpp), which its getters check
    template<typename T>
    static void samrai_type_attribute_(hid_t const dataset)
    {
        int const type_key = std::is_same_v<T, double> ? 5 : 7;
        hid_t const space  = H5Screate(H5S_SCALAR);
        check_(space, "H5Screate");
        hid_t const attr
            = H5Acreate(dataset, "Type", H5T_STD_I8BE, space, H5P_DEFAULT, H5P_DEFAULT);
        check_(attr, "H5Acreate");
        check_(H5Awrite(attr, H5T_NATIVE_INT, &type_key), "H5Awrite");
        H5Aclose(attr);
        H5Sclose(space);
    }

    template<typename T, typename Particles>
    static void put_member_(hid_t const group, std::string const& key, Particles const& particles,
                            std::size_t const offset, std::size_t const N)
    {
        hsize_t const n_values = particles.size() * N;
        hid_t const fspace     = H5Screate_simple(1, &n_values, nullptr);
        check_(fspace, "H5Screate_simple");
        hid_t const dataset = H5Dcreate(group, key.c_str(), file_type_<T>(), fspace, H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT);
        check_(dataset, "H5Dcreate " + key);
        hid_t const mspace = memspace_<T, Particles>(particles.size(), offset, N);
        check_(H5Dwrite(dataset, native_type_<T>(), mspace, H5S_ALL, H5P_DEFAULT,
                        particles.vector().data()),
               "H5Dwrite " + key);
        samrai_type_attribute_<T>(dataset);
        H5Sclose(mspace);
        H5Dclose(dataset);
        H5Sclose(fspace);
    }

    template<typename T, typename Particles>
    static void get_member_(hid_t const group, std::string const& key, Particles& particles,
                            std::size_t const offset, std::size_t const N)
    {
        hid_t const dataset = H5Dopen(group, key.c_str(), H5P_DEFAULT);
        check_(dataset, "H5Dopen " + key);
        hid_t const fspace = H5Dget_space(dataset);
        check_(fspace, "H5Dget_space " + key);
        auto const n_values = H5Sget_simple_extent_npoints(fspace);
        if (n_values < 0 or static_cast<std::size_t>(n_values) != particles.size() * N)
            throw std::runtime_error("ParticlesRestart: " + key + " has " + std::to_string(n_values)
                                     + " values, expected " + std::to_string(particles.size() * N));
        hid_t const mspace = memspace_<T, Particles>(particles.size(), offset, N);
        check_(H5Dread(dataset, native_type_<T>(), mspace, H5S_ALL, H5P_DEFAULT,
                       particles.vector().data()),
               "H5Dread " + key);
        H5Sclose(mspace);
        H5Sclose(fspace);
        H5Dclose(dataset);
    }
#else
    static int group_id_(SAMRAI::tbox::Database&)
    {
        throw std::runtime_error("ParticlesRestart: restarts need SAMRAI with HDF5");
    }
    template<typename T, typename Particles>
    static void put_member_(int, std::string const&, Particles const&, std::size_t, std::size_t)
    {
    }
    template<typename T, typename Particles>
    static void get_member_(int, std::string const&, Particles&, std::size_t, std::size_t)
    {
    }
#endif
};

} // namespace PHARE::amr

#endif /* PHARE_AMR_DATA_PARTICLES_PARTICLES_RESTART_HPP */
