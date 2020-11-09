
#ifndef PHARE_HYBRID_TAGGER_H
#define PHARE_HYBRID_TAGGER_H

#include <memory>
#include <utility>
#include <stdexcept>

#include "tagger.h"

namespace PHARE::amr
{
class HybridTaggerStrategy
{
public:
    virtual void tag()              = 0;
    virtual ~HybridTaggerStrategy() = 0;
};


HybridTaggerStrategy::~HybridTaggerStrategy() {}


class ScaledAvgHybridTaggerStrategy : public HybridTaggerStrategy
{
public:
    void tag() override;
};


void ScaledAvgHybridTaggerStrategy::tag() {}


class HybridTagger : public Tagger
{
public:
    HybridTagger(std::unique_ptr<HybridTaggerStrategy> strat)
        : strat_{std::move(strat)}
    {
    }

    void tag() override;

private:
    std::unique_ptr<HybridTaggerStrategy> strat_;
};


void HybridTagger::tag()
{
    if (strat_)
    {
        strat_->tag();
    }
    else
        throw std::runtime_error("invalid tagging strategy");
}

} // namespace PHARE::amr

#endif
