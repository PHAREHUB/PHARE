

#ifndef PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_DEF_SERIALIZING
#define PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_DEF_SERIALIZING

#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE::core
{

template<auto layout_mde, auto alloc_mde>
struct ParticlesSerializer
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<typename Src>
    void operator()(std::string const& file_name, Src const& src);

    auto open_file_to_write(std::string const& file_name) { return WriteFile{file_name}; }

    struct WriteFile
    {
        WriteFile(std::string const& file_name)
            : f{file_name, std::ios::binary}
        {
        }

        template<typename T>
        void write(T* p, std::size_t const s)
        {
            f.write(reinterpret_cast<char const*>(p), s * sizeof(T));
        }

        std::ofstream f;
    };
};



template<auto layout_mde, auto alloc_mde>
struct ParticlesDeserializer
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<typename Dst, typename Src = Dst>
    void operator()(std::string const& file_name, Dst& dst);

    template<typename Src, std::uint16_t N = 1, typename Fn>
    void readN(std::string const& file_name, Fn fn);

    auto open_file_from_start(std::string const& file_name) { return ReadFile{file_name}; }

    struct ReadFile
    {
        ReadFile(std::string const& file_name)
            : f{file_name, std::ios::binary}
        {
            f.unsetf(std::ios::skipws);
            std::streampos fileSize;
            f.seekg(0, std::ios::end);
            fileSize = f.tellg();
            f.seekg(0, std::ios::beg);
            nbytes = fileSize;
        }

        template<typename T>
        T* read(T* p, std::size_t const s)
        {
            f.read(reinterpret_cast<char*>(p), s * sizeof(T));
            return p;
        }

        auto size() const { return nbytes; }

        std::ifstream f;
        int nbytes = 0;
    };
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_DEF_SERIALIZING */
