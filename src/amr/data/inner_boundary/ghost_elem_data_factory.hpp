#ifndef PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_DATA_FACTORY_HPP
#define PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_DATA_FACTORY_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "ghost_elem_data.hpp"

#include <SAMRAI/hier/BoxGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDataFactory.h>
#include <SAMRAI/pdat/CellGeometry.h>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>

namespace PHARE::amr
{

template<std::size_t dim>
class GhostElemDataFactory : public SAMRAI::hier::PatchDataFactory
{
public:
    GhostElemDataFactory() = delete;

    GhostElemDataFactory(SAMRAI::hier::IntVector ghost, std::string const& name)
        : SAMRAI::hier::PatchDataFactory{ghost}
        , name_{name}
    {
    }

    std::shared_ptr<SAMRAI::hier::PatchDataFactory>
    cloneFactory(SAMRAI::hier::IntVector const&) final
    {
        return std::make_shared<GhostElemDataFactory>(d_ghosts, name_);
    }

    std::shared_ptr<SAMRAI::hier::PatchData> allocate(SAMRAI::hier::Patch const& patch) const final
    {
        return std::make_shared<GhostElemPatchData<dim>>(patch.getBox(), d_ghosts, name_);
    }

    std::shared_ptr<SAMRAI::hier::BoxGeometry>
    getBoxGeometry(SAMRAI::hier::Box const& box) const final
    {
        return std::make_shared<SAMRAI::pdat::CellGeometry>(box, d_ghosts);
    }

    std::size_t getSizeOfMemory(SAMRAI::hier::Box const& /*box*/) const final
    {
        throw std::runtime_error("cannot compute size from box");
    }

    bool fineBoundaryRepresentsVariable() const final { return true; }

    bool dataLivesOnPatchBorder() const final { return false; }

    bool validCopyTo(std::shared_ptr<SAMRAI::hier::PatchDataFactory> const& dst) const final
    {
        return std::dynamic_pointer_cast<GhostElemDataFactory>(dst) != nullptr;
    }

private:
    std::string const name_;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_DATA_FACTORY_HPP
