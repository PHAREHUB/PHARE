#ifndef PHARE_DEBUGOD_HPP
#define PHARE_DEBUGOD_HPP

#include "core/def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "amr/data/field/field_data.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include <SAMRAI/hier/PatchHierarchy.h>

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace PHARE::amr
{


template<typename PHARE_TYPES>
class DEBUGOD
{
public:
    static constexpr auto dimension    = PHARE_TYPES::dimension;
    static constexpr auto interp_order = PHARE_TYPES::interp_order;
    using Grid_t                       = PHARE_TYPES::Grid_t;
    using GridLayout_t                 = PHARE_TYPES::GridLayout_t;
    using FieldData_t                  = PHARE::amr::FieldData<GridLayout_t, Grid_t>;

    using Point_t = PHARE::core::Point<double, dimension>;

    struct GodValue
    {
        Point_t coords;
        std::array<int, dimension> loc_index;
        std::array<int, dimension> amr_index;
        double value;
        int rank;
        std::string patchID;

        // Add other necessary fields and methods as needed
    };
    using GodExtract = std::unordered_map<std::uint32_t, std::vector<GodValue>>;

    NO_DISCARD static DEBUGOD& INSTANCE();


    void setHierarchy(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hier)
    {
        hierarchy_ = hier;
    }

    bool isActive() const { return hierarchy_ != nullptr; }

    NO_DISCARD auto time_is(std::string name, double time) const
    {
        // take first patch of first level
        // it should be enought to get time
        // this limits to looking at times at coarser time steps for now
        auto patch = *(hierarchy_->getPatchLevel(0)->begin());
        auto pdata = getPatchData(*patch, name);
        return time == pdata->getTime();
    }

    NO_DISCARD auto inspect(std::string name, Point_t const& lower, Point_t const& upper) const
    {
        GodExtract god_values;
        for (auto ilvl = 0u; ilvl < hierarchy_->getNumberOfLevels(); ++ilvl)
        {
            auto level       = hierarchy_->getPatchLevel(ilvl);
            god_values[ilvl] = std::vector<GodValue>{};

            for (auto& patch : *level)
            {
                if (!is_local(*patch))
                    continue;

                auto extract_box = PHARE::core::Box<double, dimension>{lower, upper};
                auto patch_ghost_box
                    = phare_box_from<dimension, double>(getPatchData(*patch, name)->getGhostBox());

                auto& field    = getField(*patch, name);
                auto layout    = layoutFromPatch<GridLayout_t>(*patch);
                auto centering = GridLayout_t::centering(field.physicalQuantity());

                for (auto i = 0u; i < dimension; ++i)
                {
                    if (centering[i] == PHARE::core::QtyCentering::primal)
                    {
                        extract_box.upper[i] += 1;
                    }
                }

                auto intersected_box = patch_ghost_box * extract_box;

                if (!intersected_box)
                    continue;



                // loop on nodes
                // given the mesh_size_ on root level
                // it is easy to get the level mesh size
                // and given the lower/upper bounds of the coordinates
                // it's easy to iterate over all nodes
                // these if constexpr may be removable
                // with the FieldBox object maybe....

                GodValue gval;
                auto box = *intersected_box;

                if constexpr (dimension == 1)
                {
                    //
                }

                else if constexpr (dimension == 2)
                {
                    auto& dl     = layout.meshSize();
                    auto ixStart = static_cast<int>((box.lower[0] - layout.origin()[0]) / dl[0]);
                    auto ixEnd   = static_cast<int>((box.upper[0] - layout.origin()[0]) / dl[0]);
                    auto iyStart = static_cast<int>((box.lower[1] - layout.origin()[1]) / dl[1]);
                    auto iyEnd   = static_cast<int>((box.upper[1] - layout.origin()[1]) / dl[1]);

                    for (auto ix = ixStart; ix <= ixEnd; ++ix)
                    {
                        for (auto iy = iyStart; iy <= iyEnd; ++iy)
                        {
                            gval.coords
                                = layout.fieldNodeCoordinates(field, layout.origin(), ix, iy);
                            gval.value     = field(ix, iy);
                            gval.patchID   = to_string(patch->getGlobalId());
                            gval.rank      = get_rank(*patch);
                            gval.loc_index = {ix, iy};
                        }
                    }
                }
                else if constexpr (dimension == 3)
                {
                    // for (auto& node : intersected_box)
                    // {
                    // }
                }
                god_values[ilvl].push_back(gval);
            }
        }

        return god_values;
    }


    NO_DISCARD auto inspect(std::string name, Point_t const& coord)
    {
        return inspect(name, coord, coord);
    }



    void print(GodExtract const& god_values)
    {
        for (auto& [ilvl, values] : god_values)
        {
            std::cout << "Level " << ilvl << ":\n";
            for (auto& v : values)
            {
                auto& coords  = v.coords;
                auto& loc_idx = v.loc_index;
                auto& amr_idx = v.loc_index;
                auto& rank    = v.rank;
                auto& patchID = v.patchID;

                std::cout << "\n";
            }
        }
    }


    // void stop() { god_.release(); }

    // NO_DISCARD auto& god()
    // {
    //     if (!god_)
    //         init();
    //     return *god_;
    // }

private:
    auto get_rank(SAMRAI::hier::Patch const& patch) const
    {
        return patch.getBox().getBoxId().getOwnerRank();
    }


    bool is_local(SAMRAI::hier::Patch const& patch) const
    {
        return get_rank(patch) == PHARE::core::mpi::rank();
    }

    auto getPatchData(SAMRAI::hier::Patch const& patch, std::string const& name) const
    {
        auto db      = SAMRAI::hier::VariableDatabase::getDatabase();
        auto var_id  = db->getVariable(name);
        auto context = db->getContext("default");
        return patch.getPatchData(var_id, context);
    }

    auto& getField(SAMRAI::hier::Patch const& patch, std::string const& name) const
    {
        auto pdata            = getPatchData(patch, name);
        auto const& fielddata = std::dynamic_pointer_cast<FieldData_t>(pdata);
        return fielddata->field;
    }



    DEBUGOD() {}
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
};
}; // namespace PHARE::amr

#endif
