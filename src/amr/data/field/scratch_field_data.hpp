#ifndef PHARE_AMR_DATA_FIELD_SCRATCH_FIELD_DATA_HPP
#define PHARE_AMR_DATA_FIELD_SCRATCH_FIELD_DATA_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep

#include "amr/data/field/centered_factory.hpp"
#include "amr/data/field/field_data_factory.hpp"
#include "amr/resources_manager/field_resource.hpp"

#include <SAMRAI/hier/IntVector.h>

#include <cstddef>
#include <string>
#include <utility>

/*
 * Scratch fields for the border sums. Their data is that of View_t, with the usual ghost width,
 * but SAMRAI is told a ghost width larger by extra_ghosts. SAMRAI only links a source patch to a
 * destination patch if the source interior is within the destination ghost width (+1 for primal
 * data), so a patch two patches away, whose particles still deposit on the outermost ghost nodes,
 * is otherwise never summed in (#1575).
 */
namespace PHARE::amr
{
template<typename View_t, std::size_t extra_ghosts = 1>
class ScratchFieldUser : public View_t
{
public:
    using View_t::View_t;

    std::size_t static constexpr scratch_ghosts = extra_ghosts;
};

template<typename T>
concept is_scratch_field = requires { T::scratch_ghosts; };


// allocates and gives geometries as BaseFactory, but reports the larger ghost width to SAMRAI
template<typename BaseFactory, std::size_t extra_ghosts>
class ScratchFactory : public BaseFactory
{
public:
    template<typename... Args>
    explicit ScratchFactory(Args&&... args)
        : BaseFactory(std::forward<Args>(args)...)
    {
        this->d_ghosts
            += SAMRAI::hier::IntVector{this->d_ghosts.getDim(), static_cast<int>(extra_ghosts)};
    }
};


// the variable of Variable_t, with its factory replaced by the scratch one
template<typename Variable_t, typename BaseFactory, typename GridLayoutT, std::size_t extra_ghosts>
class ScratchVariable : public Variable_t
{
public:
    ScratchVariable(std::string const& name, auto const qty)
        : Variable_t{name, qty}
    {
        this->setPatchDataFactory(
            make_centered_factory<ScratchFactory<BaseFactory, extra_ghosts>, GridLayoutT>(
                qty, this->fineBoundaryRepresentsVariable(), this->dataLivesOnPatchBorder(), name,
                qty));
    }
};


// specialized here for fields, in scratch_tensor_field_data.hpp for tensor fields
template<typename View_t, typename Grid_t, typename GridLayoutT, typename PhysicalQuantity>
struct UserScratchFieldType;

template<typename View_t, typename Grid_t, typename GridLayoutT, typename PhysicalQuantity>
    requires(!requires { View_t::rank; })
struct UserScratchFieldType<View_t, Grid_t, GridLayoutT, PhysicalQuantity>
{
    using user_type       = UserFieldType<Grid_t, GridLayoutT>;
    using factory_type    = FieldDataFactory<GridLayoutT, Grid_t>;
    using patch_data_type = user_type::patch_data_type;
    using variable_type   = ScratchVariable<typename user_type::variable_type, factory_type,
                                            GridLayoutT, View_t::scratch_ghosts>;
};

} // namespace PHARE::amr

#endif /* PHARE_AMR_DATA_FIELD_SCRATCH_FIELD_DATA_HPP */
