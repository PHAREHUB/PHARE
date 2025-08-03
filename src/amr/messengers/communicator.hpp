#ifndef PHARE_QUANTITY_REFINER_HPP
#define PHARE_QUANTITY_REFINER_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep


#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/CoarsenSchedule.h>
#include <SAMRAI/xfer/PatchLevelBorderFillPattern.h>
#include <SAMRAI/xfer/PatchLevelInteriorFillPattern.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>


#include <map>
#include <memory>
#include <optional>

#include "core/def.hpp"


namespace PHARE
{
namespace amr
{
    struct RefinerTypes
    {
        using schedule_type  = SAMRAI::xfer::RefineSchedule;
        using algorithm_type = SAMRAI::xfer::RefineAlgorithm;
    };

    struct SynchronizerTypes
    {
        using schedule_type  = SAMRAI::xfer::CoarsenSchedule;
        using algorithm_type = SAMRAI::xfer::CoarsenAlgorithm;
    };



    template<typename ComType, std::size_t dimension>
    class Communicator
    {
    private:
        using Schedule  = typename ComType::schedule_type;
        using Algorithm = typename ComType::algorithm_type;

        std::map<int, std::map<Algorithm* const, std::shared_ptr<Schedule>>> schedules_;

        static constexpr auto is_refiner      = std::is_same_v<ComType, RefinerTypes>;
        static constexpr auto is_synchronizer = std::is_same_v<ComType, SynchronizerTypes>;


    public:
        Communicator() {}
        virtual ~Communicator() {}
        Communicator(Communicator const&) = delete;
        Communicator(Communicator&&)      = default;

        // we have an algorithm for each quantity, like Bx, By, Bz
        // even if they are to be refined/synced together.
        // the reason is that SAMRAI assumes that all Variables registered
        // to an Algorithm that are of the same type will correspond to data
        // having the same layout. This is what they do with, e.g. CellVariable
        // and CellData. This assumption ("Equivalent Classes")allow them to
        // take the geometry from the first of the quantities that share the same Variable type.
        // In PHARE all our FieldDatas have the same FieldVariable type but
        // the data layout is managed by the GridLayout and will be different.
        // For instance, although  Bx, By, Bz are all FieldData and FieldVariable
        // their data layout is different and their shape is thus different.
        std::vector<std::unique_ptr<Algorithm>> algos;

        auto& add_algorithm()
        {
            if constexpr (is_refiner)
                return algos.emplace_back(std::make_unique<SAMRAI::xfer::RefineAlgorithm>());

            if constexpr (is_synchronizer)
                return algos.emplace_back(std::make_unique<SAMRAI::xfer::CoarsenAlgorithm>(
                    SAMRAI::tbox::Dimension{dimension}));
        }


        /**
         * @brief findSchedule returns the schedule at a given level number if there is one
         * (optional).
         */
        NO_DISCARD auto& findSchedule(std::unique_ptr<Algorithm> const& algo, int levelNumber) const
        {
            if (!schedules_.count(levelNumber))
                throw std::runtime_error("no schedule for level " + std::to_string(levelNumber));

            if (!schedules_.at(levelNumber).count(algo.get()))
                throw std::runtime_error("Algorithm has not been registered with Communicator");

            return schedules_.at(levelNumber).at(algo.get());
        }



        /**
         * @brief add is used to add a refine schedule for the given level number.
         * Note that already existing schedules at this level number are overwritten.
         */
        void add(std::unique_ptr<Algorithm> const& algo, std::shared_ptr<Schedule>& schedule,
                 int levelNumber)
        {
            schedules_[levelNumber][algo.get()] = schedule;
        }

        // oevrload with  rvalue schedule
        void add(std::unique_ptr<Algorithm> const& algo, std::shared_ptr<Schedule>&& schedule,
                 int levelNumber)
        {
            add(algo, schedule, levelNumber);
        }
    };




} // namespace amr
} // namespace PHARE


#endif
