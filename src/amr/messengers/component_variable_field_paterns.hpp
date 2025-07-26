#ifndef PHARE_AMR_MESSENGERS_COMPONENT_VARIABLE_FIELD_PATTERNS_HPP
#define PHARE_AMR_MESSENGERS_COMPONENT_VARIABLE_FIELD_PATTERNS_HPP

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"

namespace PHARE::amr
{
// when registering different components to the same algorithm in SAMRAI, as we want to do for
// vecfields, we need those components not to be considered as equivalent_classes by SAMRAI.
// Without this precaution SAMRAI will assume the same geometry for all.
class XVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class YVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class ZVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

// these are mostly for refluxing where we need a different geometry for the hydro fluxes and E
class YZVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class XZVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

class XYVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
{
};

} // namespace PHARE::amr

#endif // PHARE_AMR_MESSENGERS_COMPONENT_VARIABLE_FIELD_PATTERNS_HPP
