#ifndef PHARE_DEBUGOD_HPP
#define PHARE_DEBUGOD_HPP

#include "core/def.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "core/utilities/box/box.hpp"

#include <SAMRAI/hier/PatchHierarchy.h>

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace PHARE::amr
{


template<std::size_t dim>
class DEBUGOD
{
public:
    using Point_t = PHARE::core::Point<double, dim>;

    struct GodValue
    {
        Point_t coords;
        std::array<int, dim> loc_index;
        std::array<int, dim> amr_index;
        double value;
        int rank;
        int patchID;

        // Add other necessary fields and methods as needed
    };
    using GodExtract = std::unordered_map<std::uint32_t, std::vector<GodValue>>;

    NO_DISCARD static DEBUGOD& INSTANCE();

    void init() { god_ = std::make_unique<DEBUGOD>(); }

    void setHierarchy(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hier)
    {
        hierarchy_ = hier;
    }

    NO_DISCARD auto inspect(std::string name, Point_t& lower, Point_t& upper) const
    {
        GodExtract god_values;
        for (auto ilvl = 0u; ilvl < hierarchy_->getNumberOfLevels(); ++ilvl)
        {
            auto level       = hierarchy_->getPatchLevel(ilvl);
            god_values[ilvl] = std::vector<GodValue>{};

            for (auto& patch : *level)
            {
                if (!is_local(patch))
                    continue;

                auto extract_box = PHARE::core::Box<double, dim>{lower, upper};
                auto patch_ghost_box
                    = phare_box_from<dim, double>(getPatchData(patch, name)->getGhostBox());

                auto intersected_box = *extract_box;

                if (intersected_box.isEmpty())
                    continue;

                auto& field = getField(patch, name);
                auto layout = getLayoutFromPatch(patch, field->physicalQuantity());

                // loop on nodes
                // given the mesh_size_ on root level
                // it is easy to get the level mesh size
                // and given the lower/upper bounds of the coordinates
                // it's easy to iterate over all nodes
                // these if constexpr may be removable
                // with the FieldBox object maybe....

                GodValue gval;

                if constexpr (dim == 1)
                {
                    //
                }

                else if constexpr (dim == 2)
                {
                    auto& dl     = layout.meshSize();
                    auto ixStart = intersected_box.lower(0) / dl[0] - layout.origin()[0];
                    auto ixEnd   = intersected_box.upper(0) / dl[0] - layout.origin()[0];

                    auto iyStart = intersected_box.lower(1) / dl[1] - layout.origin()[1];
                    auto iyEnd   = intersected_box.upper(1) / dl[1] - layout.origin()[1];

                    for (auto ix = ixStart; ix <= ixEnd; ++ix)
                    {
                        for (auto iy = iyStart; iy <= iyEnd; ++iy)
                        {
                            gval.coords
                                = layout.fieldNodeCoordinates(field, layout.origin(), ix, iy);
                            gval.value = field(ix, iy);
                        }
                    }
                }
                else if constexpr (dim == 3)
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


    NO_DISCARD auto inspect(std::string name, Point_t& coord)
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
    bool is_local(SAMRAI::hier::Patch const& patch) const
    {
        return patch.getBox().getBoxId().getOwnerRank() == PHARE::core::mpi::rank();
    }

    auto& getPatchData(SAMRAI::hier::Patch const& patch, std::string const& name)
    {
        auto db      = SAMRAI::hier::VariableDatabase::getDatabase();
        auto var_id  = db->getVariable(name);
        auto context = db->getContext("default");
        return patch.getPatchData(var_id, context);
    }

    auto& getField(SAMRAI::hier::Patch const& patch, std::string const& name)
    {
        auto db     = SAMRAI::hier::VariableDatabase::getDatabase();
        auto var_id = db->getVariable(name);
        // auto context   = db->getContext("default");
        // auto patchdata = patch.getPatchData(var_id, context);
        return FieldData::getPatchData(patch, var_id);
    }



    DEBUGOD() {}
    std::unique_ptr<DEBUGOD> god_;
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
};
}; // namespace PHARE::amr

#endif
