
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/CellVariable.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/pdat/CellData.h>

#include "core/data/ions/ions.h"


class WorkLoadEstimator
{
    public:
        WorkLoadEstimator() { }

        virtual void estimate(SAMRAI::hier::PatchLevel lev, double*, PHARE::core::Ions const&) = 0;
        virtual void set_strategy(std::string) = 0;
        int getID() { return id_; };

    private:
        SAMRAI::tbox::Dimension dim_{2};
        std::shared_ptr<SAMRAI::pdat::CellData> workLoad_;
        std::shared_ptr<SAMRAI::hier::CellVariable<double>> workLoadVariable_;
        std::string contextName_{"default"};
        std::shared_ptr<SAMRAI::hier::VariableContext> context_;
        SAMRAI::hier::IntVector gw0 = SAMRAI::hier::IntVector::getZero(dim_);
        int const id_{12};
};

