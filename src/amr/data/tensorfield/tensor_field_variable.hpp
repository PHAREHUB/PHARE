#ifndef PHARE_TENSORFIELD_VARIABLE_HPP
#define PHARE_TENSORFIELD_VARIABLE_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"

#include "amr/data/tensorfield/tensor_field_data_factory.hpp"

#include <SAMRAI/hier/Variable.h>

#include <utility>


namespace PHARE::amr
{

template<std::size_t rank, typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
class TensorFieldVariable : public SAMRAI::hier::Variable
{
    using tensor_t = PhysicalQuantity::template TensorType<rank>;

public:
    static constexpr std::size_t dimension    = GridLayoutT::dimension;
    static constexpr std::size_t interp_order = GridLayoutT::interp_order;

    /** \brief Construct a new variable with an unique name, and a specific PhysicalQuantity
     *
     *  TensorFieldVariable represent a data on a patch, it does not contain the data itself,
     *  after creation, one need to register it with a context : see registerVariableAndContext.
     */
    TensorFieldVariable(std::string const& name, tensor_t qty,
                        bool fineBoundaryRepresentsVariable = false)
        : SAMRAI::hier::Variable(
              name,
              std::make_shared<TensorFieldDataFactory<rank, GridLayoutT, Grid_t, PhysicalQuantity>>(
                  fineBoundaryRepresentsVariable, computeDataLivesOnPatchBorder_(qty), name, qty))
        , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
        , dataLivesOnPatchBorder_{computeDataLivesOnPatchBorder_(qty)}
    {
    }


    // The fine boundary representation boolean argument indicates which values (either coarse
    // or fine) take precedence at coarse-fine mesh boundaries during coarsen and refine
    // operations. The default is that fine data values take precedence on coarse-fine
    // interfaces.
    bool fineBoundaryRepresentsVariable() const final { return fineBoundaryRepresentsVariable_; }



    /** \brief Determines whether or not if data may lives on patch border
     *
     *  It will be true if in at least one direction, the data is primal
     */
    bool dataLivesOnPatchBorder() const final { return dataLivesOnPatchBorder_; }

private:
    bool const fineBoundaryRepresentsVariable_ = false;
    bool const dataLivesOnPatchBorder_         = false;



    bool static computeDataLivesOnPatchBorder_(tensor_t const& qty)
    {
        auto qts = PhysicalQuantity::componentsQuantities(qty);

        for (auto const& qt : qts)
        {
            auto const& centering = GridLayoutT::centering(qt);

            for (auto const& qtyCentering : centering)
            {
                if (qtyCentering == core::QtyCentering::primal)
                {
                    return true;
                }
            }
        }
        return false;
    }
};


} // namespace PHARE::amr


#endif
