#ifndef PHARE_AMR_SAMRAI_HPP
#define PHARE_AMR_SAMRAI_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/types.hpp"
#include "core/data/field/field_box.hpp"
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


// A tiled field has no single vector() of its own -- reduce it to a plain (non-tiled)
// field first, same as for diagnostics (see diagnostic::BaseModelView::field_reducer).
template<typename Grid_t>
void putFieldToRestart(auto& db, auto const& path, Grid_t const& field)
{
    if constexpr (core::is_field_tile_set_v<Grid_t>)
        putVectorToRestart(db, path, core::reduce_single(field).vector());
    else
        putVectorToRestart(db, path, field.vector());
}

template<typename Grid_t>
void getFieldFromRestart(auto& db, auto const& path, Grid_t& field)
{
    if constexpr (core::is_field_tile_set_v<Grid_t>)
    {
        typename Grid_t::grid_type grid{field.name(), field.physicalQuantity(), field.shape()};
        getVectorFromRestart(db, path, grid.vector());
        core::reduce_single(field, grid); // scatter back into tiles
    }
    else
        getVectorFromRestart(db, path, field.vector());
}


template<typename ParticleArray_t>
void putParticlesToRestart(auto& db, std::string const& name, ParticleArray_t& particles)
{
    using Packer              = core::ParticlePacker<ParticleArray_t>;
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
        putVectorToRestart(db, name + "_" + packer.keys()[part_idx++], v);
    });
}


template<typename ParticleArray_t>
void getParticlesFromRestart(auto& db, std::string const& name, ParticleArray_t& particles)
{
    using Packer              = core::ParticlePacker<ParticleArray_t>;
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
            getVectorFromRestart(db, name + "_" + Packer::keys()[part_idx++], arg);
        });
    }

    assert(particles.size() == 0);
    for (std::size_t i = 0; i < n_particles; ++i)
        particles.emplace_back(soa.copy(i));
}


} // namespace PHARE::amr

#endif /*PHARE_AMR_SAMRAI_HPP*/
