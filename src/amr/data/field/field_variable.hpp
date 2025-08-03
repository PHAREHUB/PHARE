#ifndef PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include <SAMRAI/hier/Variable.h>

#include <utility>

#include "field_data_factory.hpp"

namespace PHARE
{
namespace amr
{
    template<typename GridLayoutT, typename FieldImpl,
             typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
    /**
     * @brief The FieldVariable class
     */
    class FieldVariable : public SAMRAI::hier::Variable
    {
    public:
        static constexpr std::size_t dimension    = GridLayoutT::dimension;
        static constexpr std::size_t interp_order = GridLayoutT::interp_order;

        /** \brief Construct a new variable with an unique name, and a specific PhysicalQuantity
         *
         *  FieldVariable represent a data on a patch, it does not contain the data itself,
         *  after creation, one need to register it with a context : see registerVariableAndContext.
         *
         *
         *  Note that `fineBoundaryRepresentsVariable` is set to false so that
         *  coarse-fine interfaces are handled such that copy happens **before**
         *  refining. See https://github.com/LLNL/SAMRAI/issues/292
         */
        FieldVariable(std::string const& name, PhysicalQuantity qty,
                      bool fineBoundaryRepresentsVariable = false)
            : SAMRAI::hier::Variable(name,
                                     std::make_shared<FieldDataFactory<GridLayoutT, FieldImpl>>(
                                         fineBoundaryRepresentsVariable,
                                         computeDataLivesOnPatchBorder_(qty), name, qty))
            , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
            , dataLivesOnPatchBorder_{computeDataLivesOnPatchBorder_(qty)}
        {
        }


        // The fine boundary representation boolean argument indicates which values (either coarse
        // or fine) take precedence at coarse-fine mesh boundaries during coarsen and refine
        // operations. The default is that fine data values take precedence on coarse-fine
        // interfaces.
        bool fineBoundaryRepresentsVariable() const final
        {
            return fineBoundaryRepresentsVariable_;
        }



        /** \brief Determines whether or not if data may lives on patch border
         *
         *  It will be true if in at least one direction, the data is primal
         */
        bool dataLivesOnPatchBorder() const final { return dataLivesOnPatchBorder_; }

    private:
        bool const fineBoundaryRepresentsVariable_ = false;
        bool const dataLivesOnPatchBorder_         = false;



        bool static computeDataLivesOnPatchBorder_(PhysicalQuantity qty)
        {
            auto const& centering = GridLayoutT::centering(qty);


            for (auto const& qtyCentering : centering)
            {
                if (qtyCentering == core::QtyCentering::primal)
                {
                    return true;
                }
            }
            return false;
        }
    };
} // namespace amr


} // namespace PHARE
#endif
