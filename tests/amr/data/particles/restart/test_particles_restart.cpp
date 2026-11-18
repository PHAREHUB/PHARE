// ParticlesRestart (zero copy, HDF5 hyperslabs) against the SoA restart code it replaced, kept
//  below as legacy_put/legacy_get: files written by one must load with the other, bit for bit.
//  Files are laid out as in a SAMRAI restart (proc.0000000), for rank 0, level 0, patch 0:
//
//    /PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/
//        domainParticles_{weight,iCell,delta,v}, levelGhostParticles_*, ...
//
//  and are left in the test directory, to compare with h5ls against a real restart.

#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"

#include "amr/samrai.hpp"
#include "amr/data/particles/particles_restart.hpp"

#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/HDFDatabase.h>

#include "gtest/gtest.h"

#include <array>
#include <cassert>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
#include <stdexcept>

using namespace PHARE;


// the SoA restart code ParticlesRestart replaced, from ParticlesData::putToRestart and
//  getFromRestart before #1578
template<std::size_t dim>
void legacy_put(SAMRAI::tbox::Database& db, std::string const& name,
                core::ParticleArray<dim> const& particles)
{
    using Packer = core::ParticlePacker<dim>;

    // SAMRAI errors on writing 0 size arrays
    if (particles.size() == 0)
        return;

    Packer packer(particles);
    core::ContiguousParticles<dim> soa{particles.size()};
    packer.pack(soa);

    std::size_t part_idx = 0;
    core::apply(soa.as_tuple(), [&](auto const& arg) {
        amr::putVectorToRestart(db, name + "_" + packer.keys()[part_idx++], arg);
    });
}

template<std::size_t dim>
void legacy_get(SAMRAI::tbox::Database& db, std::string const& name,
                core::ParticleArray<dim>& particles)
{
    using Packer = core::ParticlePacker<dim>;

    std::array<bool, Packer::n_keys> keys_exist = core::generate_from(
        [&](auto const& key) { return db.keyExists(name + "_" + key); }, Packer::keys());

    bool all  = core::all(keys_exist);
    bool none = core::none(keys_exist);
    if (!(all or none))
        throw std::runtime_error("ParticlesData::getFromRestart has been given an "
                                 "invalid input file, inconsistent state detected");

    if (none) // can't read what doesn't exist
        return;

    auto n_particles = db.getArraySize(name + "_" + Packer::arbitrarySingleValueKey());
    core::ContiguousParticles<dim> soa{n_particles};

    {
        std::size_t part_idx = 0;
        core::apply(soa.as_tuple(), [&](auto& arg) {
            amr::getVectorFromRestart(db, name + "_" + Packer::keys()[part_idx++], arg);
        });
    }

    assert(particles.size() == 0);
    particles.reserve(n_particles);
    for (std::size_t i = 0; i < n_particles; ++i)
        particles.push_back(soa.copy(i));
}



// a restart file, open for writing or reading, and its patch data database
class RestartFile
{
public:
    // patch data group of rank 0, level 0, patch 0, as in a SAMRAI restart
    static constexpr std::array<char const*, 4> patch_data_path{
        "PHARE_hierarchy", "level_0000", "level_0000-patch_0000000-block_0000000",
        "protons##default"};

    RestartFile(std::string const& filename, bool const write)
    {
        if (write ? !root_->create(filename) : !root_->open(filename))
            throw std::runtime_error("cannot open " + filename);

        std::shared_ptr<SAMRAI::tbox::Database> db = root_;
        for (auto const& group : patch_data_path)
            groups_.emplace_back(db = write ? db->putDatabase(group) : db->getDatabase(group));
    }

    ~RestartFile()
    {
        groups_.clear();
        root_->close();
    }

    SAMRAI::tbox::Database& patch_data() { return *groups_.back(); }

private:
    std::shared_ptr<SAMRAI::tbox::HDFDatabase> root_
        = std::make_shared<SAMRAI::tbox::HDFDatabase>("root");
    std::vector<std::shared_ptr<SAMRAI::tbox::Database>> groups_;
};



template<typename Dim>
struct ParticlesRestartTest : public ::testing::Test
{
    static constexpr auto dim = Dim{}();
    using ParticleArray_t     = core::ParticleArray<dim>;
    using Particle_t          = ParticleArray_t::Particle_t;
    using box_t               = ParticleArray_t::box_t;

    static constexpr std::array<char const*, 3> names{"domainParticles", "levelGhostParticles",
                                                      "levelGhostParticlesOld"};

    box_t const box{core::ConstArray<int, dim>(0), core::ConstArray<int, dim>(9)};

    // every member and component distinct, so a wrong offset, stride or padding shows
    ParticleArray_t make(std::size_t const n_particles) const
    {
        ParticleArray_t particles{box};
        for (std::size_t i = 0; i < n_particles; ++i)
        {
            Particle_t p;
            p.weight = 1. / 3 + i;
            for (std::size_t d = 0; d < dim; ++d)
            {
                p.iCell[d] = static_cast<int>((i + 3 * d) % 10);
                p.delta[d] = .1 * (d + 1) + 1e-3 * i;
            }
            for (std::size_t c = 0; c < 3; ++c)
                p.v[c] = -10. * (c + 1) - 1e-2 * i;
            particles.push_back(p);
        }
        return particles;
    }

    // one array per name, the last one empty
    std::array<ParticleArray_t, 3> const particles{make(101), make(7), make(0)};

    std::string filename(std::string const& what) const
    {
        return "particles_restart_" + std::to_string(dim) + "d_" + what + ".h5";
    }

    void write(std::string const& filename, auto&& put) const
    {
        RestartFile file{filename, /*write*/ true};
        for (std::size_t i = 0; i < names.size(); ++i)
            put(file.patch_data(), names[i], particles[i]);
    }

    void read_and_check(std::string const& filename, auto&& get) const
    {
        RestartFile file{filename, /*write*/ false};
        for (std::size_t i = 0; i < names.size(); ++i)
        {
            ParticleArray_t loaded{box};
            get(file.patch_data(), names[i], loaded);
            EXPECT_EQ(loaded.size(), particles[i].size()) << names[i];
            EXPECT_TRUE(loaded == particles[i]) << names[i];
            EXPECT_EQ(loaded.nbr_particles_in(box), particles[i].size()) << names[i]; // mapped
        }
    }
};

template<std::size_t dim>
using Dim = std::integral_constant<std::size_t, dim>;

using Dims = ::testing::Types<Dim<1>, Dim<2>, Dim<3>>;
TYPED_TEST_SUITE(ParticlesRestartTest, Dims);


auto constexpr new_put = [](auto& db, auto const& name, auto const& particles) {
    amr::ParticlesRestart::put(db, name, particles);
};
auto constexpr new_get = [](auto& db, auto const& name, auto& particles) {
    amr::ParticlesRestart::get(db, name, particles);
};
auto constexpr old_put
    = [](auto& db, auto const& name, auto const& particles) { legacy_put(db, name, particles); };
auto constexpr old_get
    = [](auto& db, auto const& name, auto& particles) { legacy_get(db, name, particles); };


TYPED_TEST(ParticlesRestartTest, legacyFilesLoad)
{
    auto const filename = this->filename("legacy_write");
    this->write(filename, old_put);
    this->read_and_check(filename, new_get);
}

TYPED_TEST(ParticlesRestartTest, filesLoadWithLegacy)
{
    auto const filename = this->filename("new_write");
    this->write(filename, new_put);
    this->read_and_check(filename, old_get);
}

TYPED_TEST(ParticlesRestartTest, roundTrips)
{
    auto const filename = this->filename("round_trip");
    this->write(filename, new_put);
    this->read_and_check(filename, new_get);
}

TYPED_TEST(ParticlesRestartTest, emptyParticlesWriteNothing)
{
    auto const filename = this->filename("empty");
    this->write(filename, new_put);

    RestartFile file{filename, /*write*/ false};
    for (auto const& key : core::ParticlePacker<TestFixture::dim>::keys())
        EXPECT_FALSE(file.patch_data().keyExists(std::string{"levelGhostParticlesOld_"} + key));
}

TYPED_TEST(ParticlesRestartTest, loadingIntoNonEmptyParticlesThrows)
{
    auto const filename = this->filename("non_empty");
    this->write(filename, new_put);

    RestartFile file{filename, /*write*/ false};
    auto particles = this->make(1);
    EXPECT_THROW(new_get(file.patch_data(), "domainParticles", particles), std::runtime_error);
}

TYPED_TEST(ParticlesRestartTest, missingMembersThrow)
{
    auto const filename = this->filename("missing_members");
    {
        RestartFile file{filename, /*write*/ true};
        file.patch_data().putDoubleVector("domainParticles_weight", {1., 2.});
    }

    RestartFile file{filename, /*write*/ false};
    typename TestFixture::ParticleArray_t particles{this->box};
    EXPECT_THROW(new_get(file.patch_data(), "domainParticles", particles), std::runtime_error);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int testResult = RUN_ALL_TESTS();

    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
