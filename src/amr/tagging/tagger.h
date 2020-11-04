
#ifndef PHARE_TAGGER_H
#define PHARE_TAGGER_H

#include <memory>

namespace PHARE::amr
{
class Tagger
{
public:
    virtual void tag() = 0;
};
} // namespace PHARE::amr

#endif
