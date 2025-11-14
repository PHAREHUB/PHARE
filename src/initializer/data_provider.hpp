#ifndef DATA_PROVIDER_HPP
#define DATA_PROVIDER_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <optional>
#include <functional>

#include "cppdict/include/dict.hpp"
#include "core/def.hpp"
#include "core/utilities/span.hpp"

namespace PHARE
{
namespace initializer
{
    template<typename ReturnType, std::size_t dim>
    struct InitFunctionHelper
    {
    };

    template<>
    struct InitFunctionHelper<double, 1>
    {
        using return_type = std::shared_ptr<core::Span<double>>;
        using param_type  = std::vector<double> const&;
        using type        = std::function<return_type(param_type)>;
    };

    template<>
    struct InitFunctionHelper<double, 2>
    {
        using return_type = std::shared_ptr<core::Span<double>>;
        using param_type  = std::vector<double> const&;
        using type        = std::function<return_type(param_type, param_type)>;
    };

    template<>
    struct InitFunctionHelper<double, 3>
    {
        using return_type = std::shared_ptr<core::Span<double>>;
        using param_type  = std::vector<double> const&;
        using type        = std::function<return_type(param_type, param_type, param_type)>;
    };

    template<std::size_t dim>
    using InitFunction = typename InitFunctionHelper<double, dim>::type;


    using PHAREDict
        = cppdict::Dict<bool, int, std::vector<int>, double, std::vector<double>, std::size_t,
                        std::optional<std::size_t>, std::string, std::vector<std::string>,
                        InitFunction<1>, InitFunction<2>, InitFunction<3>>;


    class PHAREDictHandler
    {
    public:
        NO_DISCARD static PHAREDictHandler& INSTANCE();

        void init() { phareDict = std::make_unique<PHAREDict>(); }

        void stop() { phareDict.release(); }

        NO_DISCARD auto& dict()
        {
            if (!phareDict)
                init();
            return *phareDict;
        }

    private:
        PHAREDictHandler() = default;
        std::unique_ptr<PHAREDict> phareDict;
    };


    class DataProvider
    {
    public:
        /**
         * @brief read will read the data from whatever source specific DataProvider classes
         * implement. readData() must be called before passing
         */
        virtual void read()     = 0;
        virtual ~DataProvider() = default;
    };

} // namespace initializer
} // namespace PHARE

#endif // DATA_PROVIDER_HPP
