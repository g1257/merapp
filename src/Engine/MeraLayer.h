#ifndef MERALAYER_H
#define MERALAYER_H
#include "MeraFactor.h"

namespace Mera {

class MeraLayer {

    typedef MeraFactor MeraFactorType;

public:

    SizeType size() const
    {
        return data_.size();
    }

private:

    std::vector<MeraFactorType> data_;


}; //class

} //namespace

#endif // MERALAYER_H
