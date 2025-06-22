#ifndef PHARE_DEBUGOD_HPP
#define PHARE_DEBUGOD_HPP

#include "core/def.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/mpi_utils.hpp"

#include <SAMRAI/hier/PatchHierarchy.h>

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace PHARE::amr
{

struct GodValue
{
    std::vector<double> coords;
    std::vector<int> loc_index;
    std::vector<int> amr_index;
    int rank;
    int patchID;

    // Add other necessary fields and methods as needed
};

template<std::size_t dim>
class DEBUGOD
{
public:
    using Point_t    = PHARE::core::Point<double, dim>;
    using GodExtract = std::unordered_map<std::uint32_t, std::vector<GodValue>>;

    NO_DISCARD static DEBUGOD& INSTANCE();

    void init(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hier)
    {
        god_ = std::make_unique<DEBUGOD>(hier);
    }

    void setHierarchy(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hier)
    {
        hierarchy_ = hier;
    }

    NO_DISCARD inspect(std::string name, Point_t& lower, Point_t& upper)
    {
        GodExtract god_values;
        for (auto ilvl = 0u; ilvl < hierarchy_->getNumberOfLevels(); ++ilvl))
            {
                auto level       = hierarchy_->getPatchLevel(ilvl);
                god_values[ilvl] = std::vector<GodValue>{};

                for (auto& patch : *level)
                {
                    if (!is_local(patch))
                        continue;


                    auto& field      = getField(patch, name);
                    auto extract_box = Box{lower, upper};

                    for (auto& node : extract_box)
                    {
                        GodValue gval;
                        gval.coords = {node[0], node[1], node[2]};
                    }
                }
            }

        return GodExtract{};
    }


    NO_DISCARD auto inspect(std::string name, Point_t& coord)
    {
        return inspect(name, coord, coord);
    }




    void print(std::vector<GodValue> const& god_values)
    {
        for (auto& [ilvl, values] : god_values)
        {
            std::cout << "Level " << ilvl << ":\n";
            for (auto& v : values)
            {
                auto& coords  = v.coords;
                auto& loc_idx = v.loc_index;
                auto& amr_idx v.loc_index;
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
