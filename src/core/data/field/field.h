#ifndef PHARE_CORE_DATA_FIELD_FIELD_H
#define PHARE_CORE_DATA_FIELD_FIELD_H

#include <cstdint>

namespace PHARE
{
class Field
{
private:
    uint32_t size_{0};

public:
    uint32_t size() const { return size_; }
};
} // namespace PHARE
#endif
