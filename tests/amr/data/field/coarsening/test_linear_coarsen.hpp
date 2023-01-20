#ifndef PHARE_TEST_LINEAR_COARSEN_HPP
#define PHARE_TEST_LINEAR_COARSEN_HPP

#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/field/coarsening/magnetic_field_coarsener.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/field_coarsen_index_weight.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cassert>

using testing::DoubleNear;

using namespace PHARE::core;
using namespace PHARE::amr;


struct Files
{
    std::ifstream xFine, yFine, zFine;
    std::ifstream xCoarse, yCoarse, zCoarse;
    std::ifstream xCoarseAfterCoarsening, yCoarseAfterCoarsening, zCoarseAfterCoarsening;

    static Files init(std::string em, std::string dimStr)
    {
        return {std::ifstream{(em + "x_fine_original" + dimStr + "d.txt").c_str()},
                std::ifstream{(em + "y_fine_original" + dimStr + "d.txt").c_str()},
                std::ifstream{(em + "z_fine_original" + dimStr + "d.txt").c_str()},
                std::ifstream{(em + "x_coarse_original" + dimStr + "d.txt").c_str()},
                std::ifstream{(em + "y_coarse_original" + dimStr + "d.txt").c_str()},
                std::ifstream{(em + "z_coarse_original" + dimStr + "d.txt").c_str()},
                std::ifstream{(em + "x_coarse_linear_coarsed_" + dimStr + "d.txt").c_str()},
                std::ifstream{(em + "y_coarse_linear_coarsed_" + dimStr + "d.txt").c_str()},
                std::ifstream{(em + "z_coarse_linear_coarsed_" + dimStr + "d.txt").c_str()}};
    }

    bool ok() const
    {
        bool eofFlags = xFine.eof() || yFine.eof() || zFine.eof() || xCoarse.eof() || yCoarse.eof()
                        || zCoarse.eof() || xCoarseAfterCoarsening.eof()
                        || yCoarseAfterCoarsening.eof() || zCoarseAfterCoarsening.eof();

        bool badFlags = xFine.bad() || yFine.bad() || zFine.bad() || xCoarse.bad() || yCoarse.bad()
                        || zCoarse.bad() || xCoarseAfterCoarsening.bad()
                        || yCoarseAfterCoarsening.bad() || zCoarseAfterCoarsening.bad();

        return !eofFlags && !badFlags;
    }
};

template<std::size_t dim>
struct EMData
{
    using Field_t  = Field<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using FieldPtr = std::shared_ptr<Field_t>;

    std::string em_key;

    HybridQuantity::Scalar xQuantity, yQuantity, zQuantity;

    std::array<FieldPtr, 3> fineValues{};
    std::array<FieldPtr, 3> coarseValues{};
    std::array<FieldPtr, 3> expectedCoarseValues{};

    static EMData<dim> get(std::string key)
    {
        assert(key == "E" or key == "B");
        if (key == "E")
            return {"E", HybridQuantity::Scalar::Ex, HybridQuantity::Scalar::Ey,
                    HybridQuantity::Scalar::Ez};
        return {"B", HybridQuantity::Scalar::Bx, HybridQuantity::Scalar::By,
                HybridQuantity::Scalar::Bz};
    }
};



/** @brief This structure will contain data for one function to be coarsened.
 *  We get the fine and coarse data and the expected after coarse field
 *  for two centering : primal and dual
 */
template<std::size_t dimension_>
struct FieldCoarsenTestData
{
    static constexpr auto dimension  = dimension_;
    static constexpr auto interp     = 1;
    static constexpr double absError = 1.e-8;

    using GridYee_t = GridLayout<GridLayoutImplYee<dimension, interp>>;
    using Field_t   = Field<NdArrayVector<dimension>, HybridQuantity::Scalar>;

    EMData<dimension> em;

    std::array<double, dimension> meshSizeCoarse{ConstArray<double, dimension>(0.2)};
    std::array<double, dimension> meshSizeFine{ConstArray<double, dimension>(0.1)};

    std::array<int, 2> coarseIndexesX{{0, 39}};
    std::array<int, 2> fineIndexesX{18, 37};


    Point<double, dimension> coarseOrigin{
        ConstArray<double, dimension>(meshSizeCoarse[dirX] * coarseIndexesX[0])};
    Point<double, dimension> fineOrigin{
        ConstArray<double, dimension>(meshSizeFine[dirX] * fineIndexesX[0])};

    std::array<std::uint32_t, dimension> nbrCellCoarse{ConstArray<std::uint32_t, dimension>(40)};
    std::array<std::uint32_t, dimension> nbrCellFine{ConstArray<std::uint32_t, dimension>(20)};

    std::shared_ptr<GridYee_t> coarseLayout{};
    std::shared_ptr<GridYee_t> fineLayout{};


    template<typename GridLayout, typename Quantity>
    auto make_field(std::string qty_key, GridLayout const& layout, Quantity const& qty)
    {
        return std::make_shared<Field_t>(qty_key, qty, layout->allocSize(qty));
    }

    auto& init()
    {
        if constexpr (dimension == 2)
            assert(nbrCellCoarse[0] == nbrCellCoarse[1]);

        coarseLayout = std::make_shared<GridYee_t>(meshSizeCoarse, nbrCellCoarse, coarseOrigin);
        fineLayout   = std::make_shared<GridYee_t>(meshSizeFine, nbrCellFine, fineOrigin);

        em.expectedCoarseValues[0] = make_field(em.em_key + "x", coarseLayout, em.xQuantity);
        em.expectedCoarseValues[1] = make_field(em.em_key + "y", coarseLayout, em.yQuantity);
        em.expectedCoarseValues[2] = make_field(em.em_key + "z", coarseLayout, em.zQuantity);

        em.coarseValues[0] = make_field(em.em_key + "x", coarseLayout, em.xQuantity);
        em.coarseValues[1] = make_field(em.em_key + "y", coarseLayout, em.yQuantity);
        em.coarseValues[2] = make_field(em.em_key + "z", coarseLayout, em.zQuantity);

        em.fineValues[0] = make_field(em.em_key + "x", fineLayout, em.xQuantity);
        em.fineValues[1] = make_field(em.em_key + "y", fineLayout, em.yQuantity);
        em.fineValues[2] = make_field(em.em_key + "z", fineLayout, em.zQuantity);

        return *this;
    }
};

template<std::size_t dimension, typename Coarsener>
void coarsen(std::size_t idx, SAMRAI::hier::Box& fineBox, SAMRAI::hier::Box& coarseBox,
             std::array<PHARE::core::QtyCentering, dimension>& centering,
             SAMRAI::hier::IntVector& ratio, EMData<dimension>& em)
{
    Coarsener coarseIt{centering, fineBox, coarseBox, ratio};

    auto primals(ConstArray<int, dimension>());
    for (std::size_t i = 0; i < dimension; i++)
        primals[i] = centering[i] == QtyCentering::primal;

    if constexpr (dimension == 1)
    {
        for (int coarseX = 10; coarseX < 15 + primals[0]; ++coarseX)
            coarseIt(*em.fineValues[idx], *em.coarseValues[idx], Point<int, dimension>{coarseX});
    }
    else if constexpr (dimension == 2)
    {
        for (int coarseX = 10; coarseX < 15 + primals[0]; ++coarseX)
            for (int coarseY = 10; coarseY < 15 + primals[1]; ++coarseY)
                coarseIt(*em.fineValues[idx], *em.coarseValues[idx],
                         Point<int, dimension>{coarseX, coarseY});
    }
    else if constexpr (dimension == 3)
    {
        for (int coarseX = 10; coarseX < 15 + primals[0]; ++coarseX)
            for (int coarseY = 10; coarseY < 15 + primals[1]; ++coarseY)
                for (int coarseZ = 10; coarseZ < 15 + primals[2]; ++coarseZ)
                    coarseIt(*em.fineValues[idx], *em.coarseValues[idx],
                             Point<int, dimension>{coarseX, coarseY, coarseZ});
    }
};

template<std::size_t dimension, typename Coarsener>
void load(FieldCoarsenTestData<dimension>& param, Files& file_data)
{
    auto load = [](auto& value, auto const& nCells, auto& file) {
        if constexpr (dimension == 1)
        {
            for (std::uint32_t ix = 0; ix < nCells[0]; ++ix)
                file >> value(ix);
        }
        else if constexpr (dimension == 2)
        {
            for (std::uint32_t ix = 0; ix < nCells[0]; ++ix)
                for (std::uint32_t iy = 0; iy < nCells[1]; ++iy)
                    file >> value(ix, iy);
        }
        else if constexpr (dimension == 3)
        {
            for (std::uint32_t ix = 0; ix < nCells[0]; ++ix)
                for (std::uint32_t iy = 0; iy < nCells[1]; ++iy)
                    for (std::uint32_t iz = 0; iz < nCells[2]; ++iz)
                        file >> value(ix, iy, iz);
        }
    };

    auto& em = param.em;

    auto fineXNCells = param.fineLayout->allocSize(em.xQuantity);
    auto fineYNCells = param.fineLayout->allocSize(em.yQuantity);
    auto fineZNCells = param.fineLayout->allocSize(em.zQuantity);

    load(*em.fineValues[0], fineXNCells, file_data.xFine);
    load(*em.fineValues[1], fineYNCells, file_data.yFine);
    load(*em.fineValues[2], fineZNCells, file_data.zFine);

    auto coarseXNCells = param.coarseLayout->allocSize(em.xQuantity);
    auto coarseYNCells = param.coarseLayout->allocSize(em.yQuantity);
    auto coarseZNCells = param.coarseLayout->allocSize(em.zQuantity);

    load(*em.coarseValues[0], coarseXNCells, file_data.xCoarse);
    load(*em.coarseValues[1], coarseYNCells, file_data.yCoarse);
    load(*em.coarseValues[2], coarseZNCells, file_data.zCoarse);

    load(*em.expectedCoarseValues[0], coarseXNCells, file_data.xCoarseAfterCoarsening);
    load(*em.expectedCoarseValues[1], coarseYNCells, file_data.yCoarseAfterCoarsening);
    load(*em.expectedCoarseValues[2], coarseZNCells, file_data.zCoarseAfterCoarsening);

    // finally apply the coarsening
    SAMRAI::tbox::Dimension dim{dimension};

    // Here we only need the lower index so we don't have to set the upper of this boxes
    SAMRAI::hier::Box fineBoxX{dim}, fineBoxY{dim}, fineBoxZ{dim};
    SAMRAI::hier::Box coarseBoxX{dim}, coarseBoxY{dim}, coarseBoxZ{dim};

    auto centeringX = param.coarseLayout->centering(em.xQuantity);
    auto centeringY = param.coarseLayout->centering(em.yQuantity);
    auto centeringZ = param.coarseLayout->centering(em.zQuantity);

    auto setLower = [](auto& hierBox, auto& layout, auto& index, auto& centering) {
        hierBox.setLower(dirX, index - static_cast<int>(layout->nbrGhosts(centering[dirX])));

        if constexpr (dimension > 1)
            hierBox.setLower(dirY, index - static_cast<int>(layout->nbrGhosts(centering[dirY])));

        if constexpr (dimension > 2)
            hierBox.setLower(dirZ, index - static_cast<int>(layout->nbrGhosts(centering[dirZ])));
    };

    setLower(fineBoxX, param.fineLayout, param.fineIndexesX[0], centeringX);
    setLower(fineBoxY, param.fineLayout, param.fineIndexesX[0], centeringY);
    setLower(fineBoxZ, param.fineLayout, param.fineIndexesX[0], centeringZ);

    setLower(coarseBoxX, param.coarseLayout, param.coarseIndexesX[0], centeringX);
    setLower(coarseBoxY, param.coarseLayout, param.coarseIndexesX[0], centeringY);
    setLower(coarseBoxZ, param.coarseLayout, param.coarseIndexesX[0], centeringZ);

    SAMRAI::hier::IntVector ratio{dim, 2};

    coarsen<dimension, Coarsener>(0, fineBoxX, coarseBoxX, centeringX, ratio, em);
    coarsen<dimension, Coarsener>(1, fineBoxY, coarseBoxY, centeringY, ratio, em);
    coarsen<dimension, Coarsener>(2, fineBoxZ, coarseBoxZ, centeringZ, ratio, em);
}


template<std::size_t dimension>
std::vector<FieldCoarsenTestData<dimension>> createParam()
{
    std::vector<FieldCoarsenTestData<dimension>> coarsenData;

    std::string const dimStr = std::to_string(dimension);


    std::array<std::pair<std::string, Files>, 2> files{
        {{"E", Files::init("E", dimStr)}, {"B", Files::init("B", dimStr)}}};

    for (auto& [em_key, file_data] : files)
    {
        while (file_data.ok())
        {
            if (em_key == "E")
                load<dimension, DefaultFieldCoarsener<dimension>>(
                    coarsenData
                        .emplace_back(
                            FieldCoarsenTestData<dimension>{EMData<dimension>::get(em_key)})
                        .init(),
                    file_data);
            else if (em_key == "B")
            {
                load<dimension, MagneticFieldCoarsener<dimension>>(
                    coarsenData
                        .emplace_back(
                            FieldCoarsenTestData<dimension>{EMData<dimension>::get(em_key)})
                        .init(),
                    file_data);
            }
        }
    }
    return coarsenData;
}


template<typename FieldCoarsenTestDataParam>
struct FieldCoarsenOperatorTest : public ::testing::Test
{
    static constexpr auto dimension  = FieldCoarsenTestDataParam::dimension;
    static constexpr double absError = 1.e-15;
};

using FieldCoarsenOperators
    = testing::Types<FieldCoarsenTestData<1>,
                     FieldCoarsenTestData<2> /* , FieldCoarsenTestData<3>*/>; // TODO 3D
TYPED_TEST_SUITE(FieldCoarsenOperatorTest, FieldCoarsenOperators);


TYPED_TEST(FieldCoarsenOperatorTest, doTheExpectedCoarseningForEB)
{
    using FieldCoarsenTestData_t = TypeParam;
    constexpr auto dim           = FieldCoarsenTestData_t::dimension;
    constexpr auto absError      = FieldCoarsenTestData_t::absError;

    auto params = createParam<dim>();
    ASSERT_TRUE(params.size() > 0);

    auto check = [](auto& em, auto& layout, auto idx, auto& qty) {
        auto const& coarseValue         = *em.coarseValues[idx];
        auto const& expectedCoarseValue = *em.expectedCoarseValues[idx];

        auto iStartX = layout->ghostStartIndex(qty, Direction::X);
        auto iEndX   = layout->ghostEndIndex(qty, Direction::X);

        if constexpr (dim == 1)
        {
            for (auto ix = iStartX; ix <= iEndX; ++ix)
                EXPECT_THAT(coarseValue(ix), DoubleNear(expectedCoarseValue(ix), absError));
        }
        else
        {
            auto iStartY = layout->ghostStartIndex(qty, Direction::Y);
            auto iEndY   = layout->ghostEndIndex(qty, Direction::Y);

            if constexpr (dim == 2)
            {
                for (auto ix = iStartX; ix <= iEndX; ++ix)
                    for (auto iy = iStartY; iy <= iEndY; ++iy)
                        EXPECT_THAT(coarseValue(ix, iy),
                                    DoubleNear(expectedCoarseValue(ix, iy), absError));
            }
            else if constexpr (dim == 3) // TODO 3D
            {
                throw std::runtime_error("Unsupported dimension"); // uncomment/test below
                // auto iStartZ = layout->ghostStartIndex(qty, Direction::Z);
                // auto iEndZ   = layout->ghostEndIndex(qty, Direction::Z);

                // for (auto ix = iStartX; ix <= iEndX; ++ix)
                //     for (auto iy = iStartY; iy <= iEndY; ++iy)
                //         for (auto iz = iStartZ; iz <= iEndZ; ++iz)
                //             EXPECT_THAT(coarseValue(ix, iy, iz),
                //                         DoubleNear(expectedCoarseValue(ix, iy, iz), absError));
            }
        }
    };

    for (auto& param : params)
    {
        check(param.em, param.coarseLayout, 0, param.em.xQuantity);
        check(param.em, param.coarseLayout, 1, param.em.yQuantity);
        check(param.em, param.coarseLayout, 2, param.em.zQuantity);
    }
}

#endif
