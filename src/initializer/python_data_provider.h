#ifndef PHARE_PYTHON_DATA_PROVIDER_H
#define PHARE_PYTHON_DATA_PROVIDER_H

#include "data_provider.h"

#include "pragma_disable.h"

// clang-format off
DISABLE_WARNING(shadow, shadow-field-in-constructor-modified, 42)
#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/functional.h>
ENABLE_WARNING(shadow, shadow-field-in-constructor-modified, 42)
// clang-format on

namespace py = pybind11;

//#include "models/physical_state.h"
//#include <data/ions/ions.h>

namespace PHARE
{
namespace initializer
{
    extern py::scoped_interpreter guard;



    class __attribute__((visibility("hidden"))) PythonDataProvider : public DataProvider
    {
    public:
        PythonDataProvider(int argc, char const* argv)
        {
            if (argc == 2)
            {
                initModuleName_ = std::string{argv};
            }
            preparePythonPath_();
        }

        /**
         * @brief read overrides the abstract DataProvider::read method. This method basically
         * executes the user python script that fills the dictionnary.
         */
        virtual void read() override { py::eval_file(initModuleName_, scope_); }



    private:
        void preparePythonPath_() { py::eval_file("setpythonpath.py", scope_); }

        std::string initModuleName_{"job.py"};
        py::scoped_interpreter guard_{};
        py::object scope_{py::module::import("__main__").attr("__dict__")};
    };

} // namespace initializer

} // namespace PHARE

#endif // DATA_PROVIDER_H
