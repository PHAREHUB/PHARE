#ifndef PHARE_AMR_SAMRAI_HPP
#define PHARE_AMR_SAMRAI_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/types.hpp"
#include "core/vector.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"

#include <SAMRAI/tbox/RestartManager.h>
#include <SAMRAI/hier/VariableDatabase.h>
#include <SAMRAI/hier/PatchDataRestartManager.h>

#include <iostream>

namespace PHARE
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(std::string const& message, std::string const& filename, int const line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};

class SamraiLifeCycle //
{
public:
    SamraiLifeCycle(int argc = 0, char** argv = nullptr);

    ~SamraiLifeCycle();

    static void reset();

    static SAMRAI::hier::VariableDatabase* getDatabase();

    static SAMRAI::hier::PatchDataRestartManager* getPatchDataRestartManager();

    static SAMRAI::tbox::RestartManager* getRestartManager();
};


} // namespace PHARE

namespace PHARE::amr
{

template<typename T>
void getFromRestart(auto& db, auto const& path, T* data, std::size_t const size)
{
    if (size == 0)
        throw std::runtime_error("SAMRAI Restarts: vectors must be presized as expected");

    if constexpr (std::is_same_v<T, double>)
        db.getDoubleArray(path, data, size);

    else if constexpr (std::is_same_v<T, int>)
        db.getIntegerArray(path, data, size);

    else
        static_assert(core::dependent_false_v<T>,
                      "SAMRAI getFromRestart Vector not supported, add it!");
};


template<typename T, typename A, std::size_t S>
auto& getVectorFromRestart(auto& db, auto const& path, std::vector<std::array<T, S>, A>& vec)
{
    auto const size = db.getArraySize(path);
    vec.resize(size / S);
    getFromRestart(db, path, &vec[0][0], size);
    return vec;
};

template<typename T, typename A>
auto& getVectorFromRestart(auto& db, auto const& path, std::vector<T, A>& vec)
{
    auto const size = db.getArraySize(path);
    vec.resize(size);
    getFromRestart(db, path, vec.data(), size);
    return vec;
};


template<typename T>
void putToRestart(auto& db, auto const& path, T const* const data, std::size_t const size)
{
    if constexpr (std::is_same_v<T, double>)
        db.putDoubleArray(path, data, size);

    else if constexpr (std::is_same_v<T, int>)
        db.putIntegerArray(path, data, size);

    else
        static_assert(core::dependent_false_v<T>,
                      "SAMRAI putToRestart Vector not supported, add it!");
};


template<typename T, typename A, std::size_t S>
void putVectorToRestart(auto& db, auto const& path, std::vector<std::array<T, S>, A> const& vec)
{
    putToRestart(db, path, &vec[0][0], vec.size() * S);
};

template<typename T, typename A>
void putVectorToRestart(auto& db, auto const& path, std::vector<T, A> const& vec)
{
    putToRestart(db, path, vec.data(), vec.size());
};


template<typename ParticleArray_t>
void putParticlesToRestart(auto& db, std::string const& name, ParticleArray_t& particles)
{
    using Packer = core::ParticlePacker<ParticleArray_t>;
    auto constexpr static is_host_mem
        = ParticleArray_t::alloc_mode == AllocatorMode::CPU
          || ParticleArray_t::alloc_mode == AllocatorMode::GPU_UNIFIED;
    auto constexpr static dim = ParticleArray_t::dimension;

    // SAMRAI errors on writing 0 size arrays
    if (particles.size() == 0)
        return;

    if constexpr (any_in(ParticleArray_t::layout_mode, core::LayoutMode::AoSMapped))
        particles.sortMapping();

    Packer packer{particles};
    [[maybe_unused]] core::SoAParticleArray<dim> soa_;

    using enum core::LayoutMode;
    auto& soa = [&]() -> auto& {
        if constexpr (any_in(ParticleArray_t::layout_mode, AoS, AoSMapped)
                     or is_tiled(ParticleArray_t::layout_mode))
        {
            soa_.resize(particles.size());
            packer.pack(soa_);
            return soa_;
        }
        else
        {
            return particles;
        }
    }();

    std::size_t part_idx = 0;
    core::apply(soa.as_tuple(), [&](auto const& v) {
        using Vector = std::decay_t<decltype(v)>;
        using T      = typename Vector::value_type;

        auto const put_vec = [&](auto const& vec) {
            putVectorToRestart(db, name + "_" + packer.keys()[part_idx++], vec);
        };

        if constexpr (is_host_mem)
            put_vec(v);
        else
        {
            static std::vector<T> put_host_vec;
            put_host_vec.resize(v.size());
            PHARE::Vector<T>::copy(put_host_vec, v);
            put_vec(put_host_vec);
        }
    });
}


template<typename ParticleArray_t>
void getParticlesFromRestart(auto& db, std::string const& name, ParticleArray_t& particles)
{
    using Packer = core::ParticlePacker<ParticleArray_t>;
    auto constexpr static is_host_mem
        = ParticleArray_t::alloc_mode == AllocatorMode::CPU
          || ParticleArray_t::alloc_mode == AllocatorMode::GPU_UNIFIED;
    auto constexpr static dim = ParticleArray_t::dimension;

    std::array<bool, Packer::n_keys> const keys_exist = core::generate_from(
        [&](auto const& key) { return db.keyExists(name + "_" + key); }, Packer::keys());

    bool all  = core::all(keys_exist);
    bool none = core::none(keys_exist);
    if (!(all or none))
        throw std::runtime_error("getParticlesFromRestart has been given an "
                                 "invalid input file, inconsistent state detected");

    if (none) // can't read what doesn't exist
        return;

    auto n_particles = db.getArraySize(name + "_" + Packer::arbitrarySingleValueKey());
    core::SoAParticleArray<dim> soa{n_particles};

    {
        std::size_t part_idx = 0;
        core::apply(soa.as_tuple(), [&](auto& arg) {
            using Vector = std::decay_t<decltype(arg)>;
            using T      = typename Vector::value_type;

            auto& vec = [](auto& v) -> auto& {
                if constexpr (is_host_mem)
                    return v;
                else
                {
                    static std::vector<T> get_host_vec;
                    get_host_vec.clear();
                    return get_host_vec;
                }
            }(arg);

            getVectorFromRestart(db, name + "_" + Packer::keys()[part_idx++], vec);

            if constexpr (not is_host_mem)
                PHARE::Vector<T>::copy(arg, vec);
        });
    }

    assert(particles.size() == 0);
    for (std::size_t i = 0; i < n_particles; ++i)
        particles.emplace_back(soa.copy(i));

    particles.check();
}


} // namespace PHARE::amr

#endif /*PHARE_AMR_SAMRAI_HPP*/
