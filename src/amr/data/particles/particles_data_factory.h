#ifndef PHARE_PARTICLES_DATA_FACTORY_H
#define PHARE_PARTICLES_DATA_FACTORY_H

#include "particles_data.h"


#include <SAMRAI/hier/BoxGeometry.h>
#include <SAMRAI/hier/PatchDataFactory.h>
#include <SAMRAI/pdat/CellGeometry.h>



namespace PHARE
{
class ParticlesDataFactory : public SAMRAI::hier::PatchDataFactory

{
public:
    ParticlesDataFactory() = default;

    // SAMRAI interface

    // End SAMRAI interface

private:
};



} // namespace PHARE



#endif
