#ifndef PHARE_TESTS_AMR_FIELD_DATA_FIELD_DATA_TEST_PARAM_H
#define PHARE_TESTS_AMR_FIELD_DATA_FIELD_DATA_TEST_PARAM_H


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/CellDataFactory.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/pdat/NodeData.h>
#include <SAMRAI/pdat/NodeDataFactory.h>
#include <SAMRAI/pdat/NodeVariable.h>

#include "data/field/field_data.h"
#include "data/field/field_overlap.h"
#include "data/field/field_variable.h"
#include "utilities/point/point.h"

namespace PHARE
{
template<Layout layout, std::size_t dim, std::size_t interpOrder, typename FieldImpl>
struct FieldDataTestParam
{
    FieldDataTestParam(std::string const& name, HybridQuantity::Scalar quantity,
                       SAMRAI::hier::Patch& patch_0, SAMRAI::hier::Patch& patch_1)
        : field0Variable{name + std::string("_0"), true, quantity}
        , field1Variable{name + std::string("_1"), true, quantity}
        , field0Factory{field0Variable.getPatchDataFactory()}
        , field1Factory{field1Variable.getPatchDataFactory()}
        , patch0{patch_0}
        , patch1{patch_1}
        , field0Geom{field0Factory->getBoxGeometry(patch0.getBox())}
        , field1Geom{field1Factory->getBoxGeometry(patch1.getBox())}
        , field0Data{std::dynamic_pointer_cast<FieldData<layout, dim, interpOrder, FieldImpl>>(
              field0Factory->allocate(patch0))}
        , field1Data{std::dynamic_pointer_cast<FieldData<layout, dim, interpOrder, FieldImpl>>(
              field1Factory->allocate(patch1))}
    {
        auto& field0 = field0Data->field;
        auto& field1 = field1Data->field;


        auto iStart = field0Data->gridLayout.ghostStartIndex(field0, Direction::X);
        auto iEnd   = field0Data->gridLayout.ghostEndIndex(field0, Direction::X);

        for (auto ix = iStart; ix <= iEnd; ++ix)
        {
            field0(ix) = 0.0;
        }

        iStart = field1Data->gridLayout.ghostStartIndex(field1, Direction::X);
        iEnd   = field1Data->gridLayout.ghostEndIndex(field1, Direction::X);

        for (auto ix = iStart; ix <= iEnd; ++ix)
        {
            field1(ix) = 1.0;
        }
    }
    FieldVariable<layout, dim, interpOrder, FieldImpl> field0Variable;
    FieldVariable<layout, dim, interpOrder, FieldImpl> field1Variable;
    std::shared_ptr<SAMRAI::hier::PatchDataFactory> field0Factory;
    std::shared_ptr<SAMRAI::hier::PatchDataFactory> field1Factory;

    SAMRAI::hier::Patch& patch0;
    SAMRAI::hier::Patch& patch1;

    std::shared_ptr<SAMRAI::hier::BoxGeometry> field0Geom;
    std::shared_ptr<SAMRAI::hier::BoxGeometry> field1Geom;

    std::shared_ptr<FieldData<layout, dim, interpOrder, FieldImpl>> field0Data;
    std::shared_ptr<FieldData<layout, dim, interpOrder, FieldImpl>> field1Data;
};

struct Patches1D
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};


    SAMRAI::hier::Box box0{SAMRAI::hier::Index(dim, 0), SAMRAI::hier::Index(dim, 10), blockId};
    SAMRAI::hier::Box box1{SAMRAI::hier::Index(dim, 5), SAMRAI::hier::Index(dim, 20), blockId};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    double dx{0.01};
    double patch0_lo{0.};
    double patch0_hi{0.1};

    double patch1_lo{0.1};
    double patch1_hi{0.2};


    SAMRAI::hier::PatchGeometry::TwoDimBool touchesRegular{dim, false};

    std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch0Geom{
        std::make_shared<SAMRAI::geom::CartesianPatchGeometry>(SAMRAI::hier::IntVector::getOne(dim),
                                                               touchesRegular, blockId, &dx,
                                                               &patch0_lo, &patch0_hi)};

    std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch1Geom{
        std::make_shared<SAMRAI::geom::CartesianPatchGeometry>(SAMRAI::hier::IntVector::getOne(dim),
                                                               touchesRegular, blockId, &dx,
                                                               &patch1_lo, &patch1_hi)};



    SAMRAI::hier::Patch patch0{box0, patchDescriptor};
    SAMRAI::hier::Patch patch1{box1, patchDescriptor};

    Patches1D()
    {
        patch0.setPatchGeometry(patch0Geom);
        patch1.setPatchGeometry(patch1Geom);
    }
};


template<typename T>
struct AFieldData1DCenteredOnEx : public ::testing::Test
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};

    Patches1D patch1d;

    HybridQuantity::Scalar quantity{HybridQuantity::Scalar::Ex};
    std::string name{"Ex"};

    T param{name, quantity, patch1d.patch0, patch1d.patch1};


    SAMRAI::hier::IntVector ghosts{SAMRAI::hier::IntVector::getZero(this->dim)};


    std::shared_ptr<SAMRAI::pdat::CellDataFactory<double>> cell0Factory;
    std::shared_ptr<SAMRAI::pdat::CellDataFactory<double>> cell1Factory;

    std::shared_ptr<SAMRAI::pdat::CellData<double>> cell0Data;

    std::shared_ptr<SAMRAI::pdat::CellData<double>> cell1Data;
    AFieldData1DCenteredOnEx()
    {
        ghosts[0] = param.field0Data->gridLayout.nbrGhostNodes(
            param.field0Data->gridLayout.centering(quantity)[0]);
        cell0Factory = std::make_shared<SAMRAI::pdat::CellDataFactory<double>>(1, ghosts);
        cell1Factory = std::make_shared<SAMRAI::pdat::CellDataFactory<double>>(1, ghosts);

        cell0Data = std::dynamic_pointer_cast<SAMRAI::pdat::CellData<double>>(
            cell0Factory->allocate(patch1d.patch0));

        cell1Data = std::dynamic_pointer_cast<SAMRAI::pdat::CellData<double>>(
            cell1Factory->allocate(patch1d.patch1));
    }
};




template<typename T>
struct AFieldData1DCenteredOnEy : public ::testing::Test
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};

    Patches1D patch1d;

    HybridQuantity::Scalar quantity{HybridQuantity::Scalar::Ey};
    std::string name{"Ey"};

    T param{name, quantity, patch1d.patch0, patch1d.patch1};


    SAMRAI::hier::IntVector ghosts{SAMRAI::hier::IntVector::getZero(this->dim)};


    std::shared_ptr<SAMRAI::pdat::NodeDataFactory<double>> node0Factory;
    std::shared_ptr<SAMRAI::pdat::NodeDataFactory<double>> node1Factory;

    std::shared_ptr<SAMRAI::pdat::NodeData<double>> node0Data;

    std::shared_ptr<SAMRAI::pdat::NodeData<double>> node1Data;
    AFieldData1DCenteredOnEy()
    {
        ghosts[0] = param.field0Data->gridLayout.nbrGhostNodes(
            param.field0Data->gridLayout.centering(quantity)[0]);
        node0Factory = std::make_shared<SAMRAI::pdat::NodeDataFactory<double>>(1, ghosts, true);
        node1Factory = std::make_shared<SAMRAI::pdat::NodeDataFactory<double>>(1, ghosts, true);

        node0Data = std::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double>>(
            node0Factory->allocate(patch1d.patch0));

        node1Data = std::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double>>(
            node1Factory->allocate(patch1d.patch1));
    }
};

} // namespace PHARE
#endif
