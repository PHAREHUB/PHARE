#ifndef PHARE_AMR_HIERARCHY_H
#define PHARE_AMR_HIERARCHY_H

#include <algorithm>

#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/TreeLoadBalancer.h>
#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/tbox/DatabaseBox.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/MemoryDatabase.h>


#include "initializer/data_provider.h"
#include "core/utilities/point/point.h"

namespace PHARE::amr
{
/**
 * @brief The Hierarchy class is a wrapper of the SAMRAI hierarchy
 * so that users do not have to use SAMRAI types
 */
class Hierarchy : public SAMRAI::hier::PatchHierarchy
{
public:
    static auto make();
    auto const& cellWidth() const { return cellWidth_; }
    auto const& domainBox() const { return domainBox_; }
    auto const& origin() const { return origin_; }

protected:
    template<std::size_t dimension>
    Hierarchy(std::shared_ptr<SAMRAI::geom::CartesianGridGeometry>&& geo,
              std::shared_ptr<SAMRAI::tbox::MemoryDatabase>&& db,
              std::array<int, dimension> const domainBox,
              std::array<double, dimension> const origin,
              std::array<double, dimension> const cellWidth);

private:
    std::vector<double> const cellWidth_;
    std::vector<int> const domainBox_;
    std::vector<double> const origin_;
};




/**
 * @brief DimHierarchy is the concrete type of Hierarchy
 * that is built from runtime parameters in the dict.
 */
template<std::size_t _dimension>
class DimHierarchy : public Hierarchy
{
    auto shapeToBox(std::array<int, _dimension> const domainBoxShape)
    {
        auto box = domainBoxShape;
        std::for_each(std::begin(box), std::end(box), [](auto& v) { v--; });
        return box;
    }

public:
    static constexpr std::size_t dimension = _dimension;

    DimHierarchy(PHARE::initializer::PHAREDict const& dict);
};



/**
 * @brief HierarchyMaker is te functor used by makeAtRunTime to build
 * a Hierarchy.
 */
struct HierarchyMaker
{
    HierarchyMaker(PHARE::initializer::PHAREDict const& dict_);
    template<typename Dimension>
    std::shared_ptr<Hierarchy> operator()(std::size_t userDim, Dimension dimension);
    PHARE::initializer::PHAREDict const& dict;
};




//-----------------------------------------------------------------------------
//                       HierarchyMaker Definitions
//-----------------------------------------------------------------------------

inline HierarchyMaker::HierarchyMaker(PHARE::initializer::PHAREDict const& dict_)
    : dict{dict_}
{
}


template<typename Dimension>
std::shared_ptr<Hierarchy> HierarchyMaker::operator()(std::size_t userDim, Dimension dimension)
{
    if (userDim == dimension())
    {
        return std::make_shared<DimHierarchy<dimension()>>(dict);
    }
    return nullptr;
}


//-----------------------------------------------------------------------------
//                       Hierarchy Definitions
//-----------------------------------------------------------------------------


inline auto Hierarchy::make()
{
    PHARE::initializer::PHAREDict const& theDict
        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    auto dim = theDict["simulation"]["dimension"].template to<int>();
    return core::makeAtRuntime<HierarchyMaker>(dim, HierarchyMaker{theDict});
}


template<std::size_t dimension>
Hierarchy::Hierarchy(std::shared_ptr<SAMRAI::geom::CartesianGridGeometry>&& geo,
                     std::shared_ptr<SAMRAI::tbox::MemoryDatabase>&& db,
                     std::array<int, dimension> const domainBox,
                     std::array<double, dimension> const origin,
                     std::array<double, dimension> const cellWidth)
    : SAMRAI::hier::PatchHierarchy{"PHARE_hierarchy", geo, db}
    , cellWidth_(cellWidth.data(), cellWidth.data() + dimension)
    , domainBox_(domainBox.data(), domainBox.data() + dimension)
    , origin_(origin.data(), origin.data() + dimension)
{
}




//-----------------------------------------------------------------------------
//                       DimHierarchy Definitions
//-----------------------------------------------------------------------------




template<typename Type, std::size_t dimension>
void parseDimXYZType(PHARE::initializer::PHAREDict const& grid, std::string key, Type* arr)
{
    arr[0] = grid[key]["x"].template to<Type>();
    if constexpr (dimension > 1)
        arr[1] = grid[key]["y"].template to<Type>();
    if constexpr (dimension > 2)
        arr[2] = grid[key]["z"].template to<Type>();
}

template<typename Type, std::size_t dimension>
auto parseDimXYZType(PHARE::initializer::PHAREDict const& grid, std::string key)
{
    std::array<Type, dimension> arr;
    parseDimXYZType<Type, dimension>(grid, key, arr.data());
    return arr;
}

template<std::size_t dimension>
void getDomainCoords(PHARE::initializer::PHAREDict const& grid, double lower[dimension],
                     double upper[dimension])
{
    static_assert(dimension > 0 and dimension <= 3, "invalid dimension should be >0 and <=3");

    auto nbr_cells = parseDimXYZType<int, dimension>(grid, "nbr_cells");
    auto mesh_size = parseDimXYZType<double, dimension>(grid, "meshsize");
    auto origin    = parseDimXYZType<double, dimension>(grid, "origin");

    for (std::size_t i = 0; i < dimension; i++)
    {
        lower[i] = origin[i];
        upper[i] = lower[i] + nbr_cells[i] * mesh_size[i];
    }
}


template<std::size_t dimension>
auto griddingAlgorithmDatabase(PHARE::initializer::PHAREDict const& grid)
{
    static_assert(dimension > 0 and dimension <= 3, "invalid dimension should be >0 and <=3");

    auto samraiDim = SAMRAI::tbox::Dimension{dimension};
    auto db        = std::make_shared<SAMRAI::tbox::MemoryDatabase>("griddingAlgoDB");

    {
        int lowerCell[dimension], upperCell[dimension];
        std::fill_n(lowerCell, dimension, 0);
        parseDimXYZType<int, dimension>(grid, "nbr_cells", upperCell);

        upperCell[0] = grid["nbr_cells"]["x"].template to<int>() - 1;

        if constexpr (dimension >= 2)
            upperCell[1] = grid["nbr_cells"]["y"].template to<int>() - 1;

        if constexpr (dimension == 3)
            upperCell[2] = grid["nbr_cells"]["z"].template to<int>() - 1;

        std::vector<SAMRAI::tbox::DatabaseBox> dbBoxes;
        dbBoxes.push_back(SAMRAI::tbox::DatabaseBox(samraiDim, lowerCell, upperCell));
        db->putDatabaseBoxVector("domain_boxes", dbBoxes);
    }

    {
        double lowerCoord[dimension], upperCoord[dimension];
        getDomainCoords<dimension>(grid, lowerCoord, upperCoord);
        db->putDoubleArray("x_lo", lowerCoord, dimension);
        db->putDoubleArray("x_up", upperCoord, dimension);
    }

    int periodicity[dimension];
    std::fill_n(periodicity, dimension, 1); // 1==periodic, hardedcoded for all dims for now.
    db->putIntegerArray("periodic_dimension", periodicity, dimension);
    return db;
}


/*
// Required input: maximum number of levels in patch hierarchy
max_levels = 4
// Required input: vector ratio between each finer level and next coarser
ratio_to_coarser {
   level_1 = 2, 2, 2
   level_2 = 2, 2, 2
   level_3 = 4, 4, 4
}
// Optional input: int vector for largest patch size on each level.
largest_patch_size {
   level_0 = 40, 40, 40
   level_1 = 30, 30, 30
   // all finer levels will use same values as level_1...
}
// Optional input: int vector for smallest patch size on each level.
smallest_patch_size {
   level_0 = 16, 16, 16
   // all finer levels will use same values as level_0...
}
// Optional input:  buffer of one cell used on each level
proper_nesting_buffer = 1
 */
template<std::size_t dimension>
auto patchHierarchyDatabase(PHARE::initializer::PHAREDict const& amr)
{
    constexpr int ratio = 2; // Nothing else supported

    auto hierDB = std::make_shared<SAMRAI::tbox::MemoryDatabase>("HierarchyDB");

    auto maxLevelNumber = amr["max_nbr_levels"].template to<int>();
    hierDB->putInteger("max_levels", maxLevelNumber);

    std::vector<int> nesting_buffer = amr["nesting_buffer"].template to<std::vector<int>>();
    hierDB->putIntegerVector("proper_nesting_buffer", nesting_buffer);

    auto ratioToCoarserDB = hierDB->putDatabase("ratio_to_coarser");

    std::vector<int> smallestPatchSize, largestPatchSize;
    std::shared_ptr<SAMRAI::tbox::Database> smallestPatchSizeDB, largestPatchSizeDB;

    if (amr.contains("smallest_patch_size"))
    {
        smallestPatchSizeDB = hierDB->putDatabase("smallest_patch_size");
        smallestPatchSize   = amr["smallest_patch_size"].template to<std::vector<int>>();
    }

    if (amr.contains("largest_patch_size"))
    {
        largestPatchSizeDB = hierDB->putDatabase("largest_patch_size");
        largestPatchSize   = amr["largest_patch_size"].template to<std::vector<int>>();
    }

    auto addIntDimArray = [](auto& db, auto const& value, auto const& level) {
        int arr[dimension];
        std::fill_n(arr, dimension, value);
        db->putIntegerArray(level, arr, dimension);
    };

    for (auto iLevel = 0; iLevel < maxLevelNumber; ++iLevel)
    {
        std::string level{"level_" + std::to_string(iLevel)};

        if (iLevel > 0)
            addIntDimArray(ratioToCoarserDB, ratio, level);

        if (smallestPatchSizeDB)
            smallestPatchSizeDB->putIntegerVector(level, smallestPatchSize);

        if (largestPatchSizeDB)
            largestPatchSizeDB->putIntegerVector(level, largestPatchSize);
    }

    return hierDB;
}



template<std::size_t _dimension>
DimHierarchy<_dimension>::DimHierarchy(PHARE::initializer::PHAREDict const& dict)
    : Hierarchy(
        std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
            SAMRAI::tbox::Dimension{dimension}, "CartesianGridGeom",
            griddingAlgorithmDatabase<dimension>(dict["simulation"]["grid"])),
        patchHierarchyDatabase<dimension>(dict["simulation"]["AMR"]),
        shapeToBox(parseDimXYZType<int, dimension>(dict["simulation"]["grid"], "nbr_cells")),
        parseDimXYZType<double, dimension>(dict["simulation"]["grid"], "origin"),
        parseDimXYZType<double, dimension>(dict["simulation"]["grid"], "meshsize"))
{
}


} // namespace PHARE::amr
#endif
