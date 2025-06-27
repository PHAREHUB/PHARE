
#include "debug.hpp"

namespace PHARE::core
{


Debuggerino& Debuggerino::INSTANCE()
{
    static Debuggerino i;
    return i;
}


debug_scope::debug_scope(std::string const& ky)
    : key{ky}
{
    this->parent = Debuggerino::INSTANCE().stack_ptr;

    if (this->parent)
        this->parent->child = this;

    Debuggerino::INSTANCE().stack_ptr = this;
}

debug_scope::~debug_scope()
{
    Debuggerino::INSTANCE().stack_ptr = this->parent;
}


} // namespace PHARE::core
