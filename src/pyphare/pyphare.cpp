
#include "cppdict/include/dict.hpp"
#include "python_data_provider.h"

#if !defined(PHARE_PYBINDED)
#define PHARE_PYBINDED
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#endif

#include "add.h"
#include "diagnostic.h"

PYBIND11_MODULE(pyphare, m)
{
    PHARE::pybind::add(m);
    PHARE::pybind::diagnostic(m);
}
