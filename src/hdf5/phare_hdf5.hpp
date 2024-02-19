#ifndef PHARE_HDF5_PHARE_HDF5_HPP
#define PHARE_HDF5_PHARE_HDF5_HPP


#if PHARE_HAS_HIGHFIVE

#include "highfive/H5Version.hpp"

#ifndef HIGHFIVE_VERSION_STRING
#pragma message("Highfive must be at least version 2.7.1")
#error // highfive is too old
#endif

#define _PHARE_WITH_HIGHFIVE(...) __VA_ARGS__
#else
#define _PHARE_WITH_HIGHFIVE(...)
#endif



#endif // PHARE_HDF5_PHARE_HDF5_HPP
