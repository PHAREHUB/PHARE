
#ifndef PHARE_TAGGER_H
#define PHARE_TAGGER_H

#include <memory>

namespace PHARE::amr
{
class Tagger
{
public:
    virtual void tag() = 0;
    virtual ~Tagger()  = 0;
};

Tagger::~Tagger() {}

} // namespace PHARE::amr

#endif
