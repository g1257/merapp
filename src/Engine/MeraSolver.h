#ifndef MERASOLVER_H
#define MERASOLVER_H
#include <iostream>
#include "MeraStep.h"


namespace Mera {

class MeraSolver {

    typedef MeraStep MeraStepType;

public:

    void computeGroundState()
    {
        SizeType qiter = 1;
        SizeType qlayer = 1;

        for (SizeType iter=0; iter<qiter; ++iter) {
            for (SizeType layer=0; layer<qlayer; ++layer) {
                step_.optimize(iter,layer);
            }
        }


    }

private:

    MeraStepType step_;
}; //class

} //namespace

#endif
