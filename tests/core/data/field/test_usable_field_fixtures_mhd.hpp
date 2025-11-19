#ifndef PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_MHD_HPP
#define PHARE_TEST_CORE_DATA_TEST_FIELD_FIXTURES_MHD_HPP


#include "core/mhd/mhd_quantities.hpp"
#include "core/data/field/field.hpp"
#include "core/data/grid/grid.hpp"

namespace PHARE::core
{

template<std::size_t dim>
class UsableFieldMHD : public Field<dim, MHDQuantity::Scalar, double>
{
public:
    auto static constexpr dimension = dim;
    using Super                     = Field<dim, MHDQuantity::Scalar, double>;
    using Grid_t                    = Grid<NdArrayVector<dim>, MHDQuantity::Scalar>;
    using scalar_t                  = MHDQuantity::Scalar;

    template<typename GridLayout>
    UsableFieldMHD(std::string const& name, GridLayout const& layout, scalar_t qty)
        : Super{name, qty}
        , xyz{name, layout, qty}
    {
        super().setBuffer(&xyz);
    }
    void set_on(Super& field) { field.setBuffer(&xyz); }

    Super& super() { return *this; }
    Super const& super() const { return *this; }

protected:
    Grid_t xyz;
};

} // namespace PHARE::core

#endif
