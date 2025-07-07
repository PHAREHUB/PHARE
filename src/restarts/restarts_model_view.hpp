#ifndef RESTART_MODEL_VIEW_HPP
#define RESTART_MODEL_VIEW_HPP


#include "core/def.hpp"

#include <string>

namespace PHARE::restarts
{
class IModelView
{
public:
    inline virtual ~IModelView();
};
IModelView::~IModelView() {}



template<typename Hierarchy, typename Model>
class ModelView : public IModelView
{
public:
    ModelView(Hierarchy const& hierarchy, Model const& model)
        : model_{model}
        , hierarchy_{hierarchy}
    {
    }



    auto writeRestartFile(std::string const& path) const
    {
        return hierarchy_.writeRestartFile(path);
    }

    NO_DISCARD auto static restartFilePathForTime(std::string path, double timestamp)
    {
        return Hierarchy::restartFilePathForTime(path, timestamp);
    }


    NO_DISCARD auto patch_data_ids() const { return model_.patch_data_ids(); }


    ModelView(ModelView const&)            = delete;
    ModelView(ModelView&&)                 = delete;
    ModelView& operator=(ModelView const&) = delete;
    ModelView& operator=(ModelView&&)      = delete;

protected:
    Model const& model_;
    Hierarchy const& hierarchy_;
};



} // namespace PHARE::restarts



#endif // RESTART_MODEL_VIEW_HPP
