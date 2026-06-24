
#ifndef PHARE_VTK_HDF5_VERSION
#define PHARE_VTK_HDF5_VERSION 1
#endif


#if PHARE_VTK_HDF5_VERSION == 0
#include "diagnostic/detail/vtk_types/spec_type_writer.hpp" // IWYU pragma: keep
#elif PHARE_VTK_HDF5_VERSION == 1
#include "diagnostic/detail/vtk_types/hax_type_writer.hpp" // IWYU pragma: keep
#else
#error // no impl!
#endif
