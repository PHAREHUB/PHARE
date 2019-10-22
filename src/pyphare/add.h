#ifndef PHARE_PYPHARE_ADD_H
#define PHARE_PYPHARE_ADD_H

using PHARE::initializer::ScalarFunction;

namespace PHARE
{
namespace pybind
{
    template<std::size_t dim, typename T,
             typename = std::enable_if_t<cppdict::is_in<T, int, double, std::string,
                                                        ScalarFunction<1> /*, VectorFunction<1>*/>>>
    void add(std::string path, T&& value)
    {
        if constexpr (dim == 1)
        {
            std::cout << path + "\n";
            cppdict::add(path, std::forward<T>(value),
                         PHARE::initializer::PHAREDictHandler::INSTANCE().dict<dim>());
        }
        /*if constexpr (dim == 2)
        {
            std::cout << path + "\n";
            cppdict::add(path, std::forward<T>(value),
                         PHARE::initializer::PHAREDictHandler::INSTANCE().dict<dim>());
        }
        if constexpr (dim == 3)
        {
            std::cout << path + "\n";
            cppdict::add(path, std::forward<T>(value),
                         PHARE::initializer::PHAREDictHandler::INSTANCE().dict<dim>());
        }*/
    }

    template<typename PyBindModule>
    void add(PyBindModule& m)
    {
        m.def("add", add<1, int, void>, "add");
        // m.def("add", add<2, int, void>, "add");
        // m.def("add", add<3, int, void>, "add");

        m.def("add", add<1, double, void>, "add");
        // m.def("add", add<2, double, void>, "add");
        // m.def("add", add<3, double, void>, "add");


        m.def("add", add<1, std::string, void>, "add");
        // m.def("add", add<2, std::string, void>, "add");
        // m.def("add", add<3, std::string, void>, "add");

        m.def("addScalarFunction1D", add<1, ScalarFunction<1>, void>, "addScalarFunction1D");
        // m.def("add", add<2, ScalarFunction<2>, void>, "add");
        // m.def("add", add<3, ScalarFunction<3>, void>, "add");
        //
    }
} /* namespace pybind*/
} /* namespace PHARE*/

#endif /* PHARE_PYPHARE_ADD_H */
