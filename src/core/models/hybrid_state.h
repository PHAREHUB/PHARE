#ifndef HYBRID_HYBRID_STATE_H
#define HYBRID_HYBRID_STATE_H


#include "data/ions/ion_initializer.h"
#include "data_provider.h"
#include "models/physical_state.h"

#include <cstddef>
#include <utility>

namespace PHARE
{
namespace core
{
    template<typename IonsInitializer>
    class HybridStateInitializer : public PhysicalStateInitializer
    {
    public:
        IonsInitializer ionInitializer;
    };

    /**
     * @brief The HybridState class is a concrete implementation of a IPhysicalState.
     * It holds an Electromag and Ion object manipulated by Hybrid concrete type of ISolver
     */
    template<typename Electromag, typename Ions, typename IonsInitializer>
    class HybridState : public IPhysicalState
    {
    public:
        static constexpr auto dimension = Ions::dimension;

        [[deprecated]] HybridState(IonsInitializer ionsInitializer)
            : ions{std::move(ionsInitializer)}
        {
        }


        HybridState(PHARE::initializer::PHAREDict<dimension> dict)
            : ions{dict["ions"]}
        {
        }

        Electromag electromag{"EM"};
        Ions ions;




        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------

        bool isUsable() const { return electromag.isUsable() and ions.isUsable(); }



        bool isSettable() const { return electromag.isSettable() and ions.isSettable(); }


        auto getCompileTimeResourcesUserList() const
        {
            return std::forward_as_tuple(electromag, ions);
        }

        auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(electromag, ions); }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------
    };


} // namespace core
} // namespace PHARE


#endif // PHARE_HYBRID_STATE_H
