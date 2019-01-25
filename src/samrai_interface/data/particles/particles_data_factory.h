#ifndef PHARE_PARTICLES_DATA_FACTORY_H
#define PHARE_PARTICLES_DATA_FACTORY_H

#include "particles_data.h"


#include <SAMRAI/hier/BoxGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDataFactory.h>
#include <SAMRAI/pdat/CellGeometry.h>



namespace PHARE
{
namespace amr_interface
{
    template<std::size_t dim>
    /**
     * @brief The ParticlesDataFactory class
     */
    class ParticlesDataFactory : public SAMRAI::hier::PatchDataFactory

    {
    public:
        ParticlesDataFactory() = delete;

        // SAMRAI interface

        ParticlesDataFactory(SAMRAI::hier::IntVector ghost, bool fineBoundaryRepresentsVariable)
            : SAMRAI::hier::PatchDataFactory{ghost}
            , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
        {
        }

        virtual std::shared_ptr<SAMRAI::hier::PatchDataFactory>
        cloneFactory(SAMRAI::hier::IntVector const&) final
        {
            return std::make_shared<ParticlesDataFactory>(d_ghosts,
                                                          fineBoundaryRepresentsVariable_);
        }

        virtual std::shared_ptr<SAMRAI::hier::PatchData>
        allocate(const SAMRAI::hier::Patch& patch) const final
        {
            return std::make_shared<ParticlesData<dim>>(patch.getBox(), d_ghosts);
        }

        virtual std::shared_ptr<SAMRAI::hier::BoxGeometry>
        getBoxGeometry(const SAMRAI::hier::Box& box) const final
        {
            return std::make_shared<SAMRAI::pdat::CellGeometry>(box, d_ghosts);
        }

        virtual size_t getSizeOfMemory(const SAMRAI::hier::Box& box) const final
        {
            throw std::runtime_error("cannot compute size from box");
        }

        virtual bool fineBoundaryRepresentsVariable() const final
        {
            return fineBoundaryRepresentsVariable_;
        }

        virtual bool dataLivesOnPatchBorder() const final { return true; }

        virtual bool validCopyTo(std::shared_ptr<SAMRAI::hier::PatchDataFactory> const&
                                     destinationPatchDataFactory) const final
        {
            auto casted
                = std::dynamic_pointer_cast<ParticlesDataFactory>(destinationPatchDataFactory);
            return casted != nullptr;
        }

        // End SAMRAI interface
    private:
        bool fineBoundaryRepresentsVariable_;


    private:
    };
} // namespace amr_interface


} // namespace PHARE



#endif
