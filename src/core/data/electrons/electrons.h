#ifndef PHARE_ELECTRONS_H
#define PHARE_ELECTRONS_H

#include "initializer/data_provider.h"

namespace PHARE::core {

template<typename Ions, typename Electromag>
class Electrons {

    using VecField = typename Ions::vecfield_type;
    using Field = typename Ions::field_type;

public:
    Electrons(PHARE::initializer::PHAREDict dict, Ions &ions, Electromag& electromag, VecField& J)
        :ions_{ions}
        ,electromag_{electromag}
        , J_{J}
    {}

    void update(){}

    bool isUsable() const
    {
        return Ne_ != nullptr && ions_.isUsable() && electromag_.isUsable() && J_.isUsable();
    }


    Field& density() {return *Ne_;}


private:
    Field* Ne_;
    Ions& ions_;
    Electromag& electromag_;
    VecField& J_;
};

}


#endif
