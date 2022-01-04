#ifndef PHARE_GRIDLAYOUT_UTILS_H
#define PHARE_GRIDLAYOUT_UTILS_H

#include <tuple>
#include <stdexcept>

namespace PHARE::core
{
template<typename GridLayout>
class LayoutHolder
{
protected:
    GridLayout* layout_{nullptr};

public:
    void setLayout(GridLayout* ptr) { layout_ = ptr; }

    bool hasLayout() const { return layout_ != nullptr; }
};




template<typename GridLayout, typename... GridLayoutSettable>
class SetLayout
{
public:
    SetLayout(GridLayout* ptr, GridLayoutSettable&... settables)
        : settables_{settables...}
    {
        std::apply([ptr](auto&... settable) { (settable.setLayout(ptr), ...); }, settables_);
    }

    ~SetLayout()
    {
        std::apply([](auto&... settable) { (settable.setLayout(nullptr), ...); }, settables_);
    }


private:
    std::tuple<GridLayoutSettable&...> settables_;
};
} // namespace PHARE::core



#endif
