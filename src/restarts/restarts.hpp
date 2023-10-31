#ifndef RESTART_RESTARTS_HPP
#define RESTART_RESTARTS_HPP

#if !defined(PHARE_HAS_HIGHFIVE)
#error // PHARE_HAS_HIGHFIVE expected to be defined as bool
#endif

#include <memory>
#include "cppdict/include/dict.hpp"

#include "hdf5/phare_hdf5.hpp"
#include "restarts_manager.hpp"


#if PHARE_HAS_HIGHFIVE

#include "restarts_model_view.hpp"
#include "restarts/detail/h5writer.hpp"

#endif

namespace PHARE::restarts
{
struct NullOpRestartsManager : public IRestartsManager
{
    void dump(double /*timeStamp*/, double /*timeStep*/) override
    {
        throw std::runtime_error("NOOP");
    }
};

struct RestartsManagerResolver
{
    template<typename Hierarchy, typename Model>
    [[nodiscard]] static std::unique_ptr<IRestartsManager>
    make_unique(Hierarchy& hier, Model& model, initializer::PHAREDict const& dict)
    {
#if PHARE_HAS_HIGHFIVE
        using ModelView_t = ModelView<Hierarchy, Model>;
        using Writer_t    = h5::Writer<ModelView_t>;
        return RestartsManager<Writer_t>::make_unique(hier, model, dict);
#else
        return std::make_unique<NullOpRestartsManager>();
#endif
    }
};

} // namespace PHARE::restarts

#endif // RESTART_RESTARTS_HPP
