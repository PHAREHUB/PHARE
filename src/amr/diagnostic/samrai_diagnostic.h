#ifndef PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H
#define PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H

#include "phare/diagnostic_manager.h"

namespace PHARE
{
namespace diagnostic
{
    struct Constants
    {
        static std::unordered_map<std::string, PHARE::core::Component> const COMPONENTS;
        static std::unordered_map<std::string, PHARE::core::Direction> const DIRECTIONS;
    };

    class KeyVector
    {
    public:
        using T = double;

        KeyVector(char const* key)
            : k(key)
        {
        }

        KeyVector(std::string&& key)
            : k(key)
        {
        }

        std::string k;
        std::vector<T> v;
    };
    template<typename MODEL>
    class SamraiPatch : public Constants
    {
    public:
        using T     = double;
        using C     = PHARE::core::Component;
        using D     = PHARE::core::Direction;
        using G     = typename MODEL::gridLayout_type;
        using R     = typename MODEL::resources_manager_type;
        using GUARD = amr::ResourcesGuard<R, ElectromagT, IonsT>;

        SamraiPatch(SAMRAI::hier::Patch& p, MODEL& m)
            : p_(p)
            , m_(m)
            , guard(m_.resourcesManager->setOnPatch(p, m_.state.electromag, m_.state.ions))
            , grid(PHARE::amr::layoutFromPatch<typename MODEL::gridLayout_type>(p))
        {
        }

        template<typename F>
        std::vector<KeyVector>& getBorE(F& f, std::vector<KeyVector>& kvs)
        {
            for (auto& kv : kvs)
            {
                getField(f, kv.v, COMPONENTS.at(kv.k), DIRECTIONS.at("x"));
            }
            return kvs;
        }

    private:
        template<typename F, typename V>
        void getField(F& f, V& v, C c, D d)
        {
            auto& field = f.getComponent(c);
            auto iStart = grid.physicalStartIndex(field, d);
            auto iEnd   = grid.physicalEndIndex(field, d);
            v.resize(iEnd - iStart + 1);
            size_t field_idx = 0;
            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                v[field_idx++] = field(ix);
            }
        }

        SAMRAI::hier::Patch& p_;
        MODEL& m_;
        GUARD guard;
        G grid;
    };
} // namespace diagnostic

std::unordered_map<std::string, PHARE::core::Component> const diagnostic::Constants::COMPONENTS
    = {{"x", PHARE::core::Component::X},
       {"y", PHARE::core::Component::Y},
       {"z", PHARE::core::Component::Z}};
std::unordered_map<std::string, PHARE::core::Direction> const diagnostic::Constants::DIRECTIONS
    = {{"x", PHARE::core::Direction::X},
       {"y", PHARE::core::Direction::Y},
       {"z", PHARE::core::Direction::Z}};

namespace amr
{
    template<typename MODEL>
    class SamraiDiagnostic : public ModelDiagnostic<MODEL>, public diagnostic::Constants
    {
    public:
        using T  = double;
        using H  = SAMRAI::hier::PatchHierarchy;
        using P  = SAMRAI::hier::Patch;
        using PL = std::shared_ptr<SAMRAI::hier::PatchLevel>;

        template<typename D, diagnostic::Level L>
        void dump(std::vector<std::shared_ptr<Diagnostic>> const&);

        auto& getLevelAt(size_t index) {}
        auto& getPatchAt(/*Level index*/) {}
        auto& getPatches(/*Level index*/) {}

    protected:
        SamraiDiagnostic(H& h, MODEL& m)
            : ModelDiagnostic<MODEL>(m)
            , hierarchy_(h)
        {
        }

        H& hierarchy_;

    private:
        SamraiDiagnostic(const SamraiDiagnostic&)             = delete;
        SamraiDiagnostic(const SamraiDiagnostic&&)            = delete;
        SamraiDiagnostic& operator&(const SamraiDiagnostic&)  = delete;
        SamraiDiagnostic& operator&(const SamraiDiagnostic&&) = delete;
    };
} /*namespace amr*/
} // namespace PHARE

#if defined(PHARE_WITH_HIGHFIVE)
#include "diagnostic/detail/samrai_highfive.h"
#endif
namespace PHARE
{
template<typename MODEL>
template<typename D, diagnostic::Level L>
void amr::SamraiDiagnostic<MODEL>::dump(std::vector<std::shared_ptr<Diagnostic>> const& diagnostics)
{
#if defined(PHARE_WITH_HIGHFIVE)
    if constexpr (std::is_same<D, hi5::Diagnostic>::value)
    {
        static_cast<SamraiHighFiveDiagnostic<MODEL>>(this)->template dump<L>(diagnostics);
    }
#endif
    throw std::runtime_error("Internal error: Unhandled Diagnostic Implementation");
}
} /*namespace PHARE*/

#endif /*PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H*/
