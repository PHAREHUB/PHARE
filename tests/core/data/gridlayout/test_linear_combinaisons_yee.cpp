
#include "test_linear_combinaisons_yee.h"

using GridLayout1DO1 = GridLayout<GridLayoutImplYee<1, 1>>;
using GridLayout1DO2 = GridLayout<GridLayoutImplYee<1, 2>>;
using GridLayout1DO3 = GridLayout<GridLayoutImplYee<1, 3>>;

using GridLayout2DO1 = GridLayout<GridLayoutImplYee<2, 1>>;
using GridLayout2DO2 = GridLayout<GridLayoutImplYee<2, 2>>;
using GridLayout2DO3 = GridLayout<GridLayoutImplYee<2, 3>>;

using GridLayout3DO1 = GridLayout<GridLayoutImplYee<3, 1>>;
using GridLayout3DO2 = GridLayout<GridLayoutImplYee<3, 2>>;
using GridLayout3DO3 = GridLayout<GridLayoutImplYee<3, 3>>;




TEST(MomentToEx, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_momentToEx.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::momentsToEx();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::momentsToEx();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::momentsToEx();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::momentsToEx();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::momentsToEx();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::momentsToEx();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::momentsToEx();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::momentsToEx();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::momentsToEx();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}



TEST(MomentToEy, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_momentToEy.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::momentsToEy();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::momentsToEy();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::momentsToEy();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::momentsToEy();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::momentsToEy();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::momentsToEy();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::momentsToEy();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::momentsToEy();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::momentsToEy();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}



TEST(MomentToEz, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_momentToEz.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::momentsToEz();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::momentsToEz();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::momentsToEz();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::momentsToEz();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::momentsToEz();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::momentsToEz();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::momentsToEz();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::momentsToEz();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::momentsToEz();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}



TEST(ExToMoment, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_ExToMoment.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::ExToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::ExToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::ExToMoments();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::ExToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::ExToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::ExToMoments();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::ExToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::ExToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::ExToMoments();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(EyToMoment, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_EyToMoment.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::EyToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::EyToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::EyToMoments();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::EyToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::EyToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::EyToMoments();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::EyToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::EyToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::EyToMoments();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(EzToMoment, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_EzToMoment.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::EzToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::EzToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::EzToMoments();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::EzToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::EzToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::EzToMoments();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::EzToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::EzToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::EzToMoments();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(JxToMoment, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_JxToMoment.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::JxToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::JxToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::JxToMoments();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::JxToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::JxToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::JxToMoments();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::JxToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::JxToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::JxToMoments();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(JyToMoment, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_JyToMoment.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::JyToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::JyToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::JyToMoments();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::JyToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::JyToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::JyToMoments();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::JyToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::JyToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::JyToMoments();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(JzToMoment, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_JzToMoment.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::JzToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::JzToMoments();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::JzToMoments();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::JzToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::JzToMoments();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::JzToMoments();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::JzToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::JzToMoments();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::JzToMoments();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(ByToEx, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_ByToEx.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::ByToEx();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::ByToEx();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::ByToEx();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::ByToEx();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::ByToEx();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::ByToEx();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::ByToEx();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::ByToEx();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::ByToEx();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(BzToEx, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_BzToEx.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::BzToEx();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::BzToEx();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::BzToEx();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::BzToEx();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::BzToEx();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::BzToEx();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::BzToEx();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::BzToEx();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::BzToEx();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(BxToEy, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_BxToEy.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::BxToEy();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::BxToEy();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::BxToEy();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::BxToEy();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::BxToEy();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::BxToEy();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::BxToEy();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::BxToEy();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::BxToEy();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(BzToEy, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_BzToEy.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::BzToEy();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::BzToEy();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::BzToEy();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::BzToEy();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::BzToEy();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::BzToEy();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::BzToEy();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::BzToEy();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::BzToEy();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}




TEST(BxToEz, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_BxToEz.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::BxToEz();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::BxToEz();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::BxToEz();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::BxToEz();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::BxToEz();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::BxToEz();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::BxToEz();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::BxToEz();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::BxToEz();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}



TEST(ByToEz, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_ByToEz.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
                EXPECT_DOUBLE_EQ(combi.coef, actual[iPoint].coef);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::ByToEz();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::ByToEz();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::ByToEz();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::ByToEz();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::ByToEz();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::ByToEz();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::ByToEz();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::ByToEz();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::ByToEz();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}



TEST(JxToEx, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_JxToEx.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::JxToEx();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::JxToEx();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::JxToEx();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::JxToEx();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::JxToEx();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::JxToEx();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::JxToEx();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::JxToEx();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::JxToEx();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}



TEST(JyToEy, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_JyToEy.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::JyToEy();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::JyToEy();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::JyToEy();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::JyToEy();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::JyToEy();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::JyToEy();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::JyToEy();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::JyToEy();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::JyToEy();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}



TEST(JzToEz, combinaisonOk)
{
    auto expectedCombinaisons = readFile("linear_coefs_yee_JzToEz.txt");
    for (auto const& combi : expectedCombinaisons)
    {
        auto dim         = combi.dimension;
        auto interpOrder = combi.interpOrder;

        auto f1D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
            }
        };
        auto f2D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
            }
        };
        auto f3D = [&combi](auto const& actual) {
            EXPECT_EQ(combi.ix.size(), actual.size());
            EXPECT_EQ(combi.iy.size(), actual.size());
            EXPECT_EQ(combi.iz.size(), actual.size());
            for (std::size_t iPoint = 0; iPoint < actual.size(); ++iPoint)
            {
                EXPECT_EQ(combi.ix[iPoint], actual[iPoint].indexes[0]);
                EXPECT_EQ(combi.iy[iPoint], actual[iPoint].indexes[1]);
                EXPECT_EQ(combi.iz[iPoint], actual[iPoint].indexes[2]);
            }
        };

        switch (dim)
        {
            case 1:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout1DO1::JzToEz();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout1DO2::JzToEz();
                    f1D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout1DO3::JzToEz();
                    f1D(actualCombinaison);
                }
                break;
            case 2:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout2DO1::JzToEz();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout2DO2::JzToEz();
                    f2D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout2DO3::JzToEz();
                    f2D(actualCombinaison);
                }
                break;
            case 3:
                if (interpOrder == 1)
                {
                    auto actualCombinaison = GridLayout3DO1::JzToEz();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 2)
                {
                    auto actualCombinaison = GridLayout3DO2::JzToEz();
                    f3D(actualCombinaison);
                }
                else if (interpOrder == 3)
                {
                    auto actualCombinaison = GridLayout3DO3::JzToEz();
                    f3D(actualCombinaison);
                }
                break;
        }
    }
}
