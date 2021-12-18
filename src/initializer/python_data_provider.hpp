
// this file includes python which has conflicting #defines with SAMRAI
//  so include it last

#ifndef PHARE_PYTHON_DATA_PROVIDER_H
#define PHARE_PYTHON_DATA_PROVIDER_H

#include "initializer/data_provider.hpp"

#include "pragma_disable.hpp"

// clang-format off
DISABLE_WARNING(shadow, shadow-field-in-constructor-modified, 42)

#undef HAVE_SYS_TIMES_H // included in python again, possibly with different value
#undef HAVE_UNISTD_H

#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/functional.h>
ENABLE_WARNING(shadow, shadow-field-in-constructor-modified, 42)
// clang-format on

namespace py = pybind11;

namespace PHARE
{
namespace initializer
{
    class __attribute__((visibility("hidden"))) PythonDataProvider : public DataProvider
    {
    public:
        PythonDataProvider() {}
        PythonDataProvider(std::string moduleName)
            : moduleName_{moduleName}
        {
        }

        /**
         * @brief read overrides the abstract DataProvider::read method. This method basically
         * executes the user python script that fills the dictionnary.
         */
        void read() override
        {
            auto module = py::module::import(initModuleName_.c_str());
            module.attr("get_user_inputs")(moduleName_);
        }



    private:
        std::string moduleName_{"job"};
        std::string initModuleName_{"pyphare.pharein.init"};
        py::scoped_interpreter guard_;
    };

} // namespace initializer

} // namespace PHARE

#endif // DATA_PROVIDER_H
