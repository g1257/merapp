#ifndef MERALAYER_H
#define MERALAYER_H
#include <cassert>
#include "Vector.h"

namespace Mera {

class MeraLayer {

    typedef int TensorType;

    enum TensorTypeEnum {TENSOR_TYPE_W,TENSOR_TYPE_U};

public:

    typedef std::pair<SizeType,SizeType> PairSizeType;

    MeraLayer()
        : vec_(u_.size())
    {
        SizeType n = vec_.size();
        SizeType nOverTwo = n/2;

        for (SizeType i=0; i<nOverTwo; ++i) {
            vec_[2*i] = PairSizeType(i,TENSOR_TYPE_W);
            vec_[2*i+1] = PairSizeType(i,TENSOR_TYPE_U);
        }
    }

    SizeType size() const
    {
        return vec_.size();
    }

    const PairSizeType& tensorToOptimize(SizeType i) const
    {
        //std::assert(i < vec_.size());
        return vec_[i];
    }

private:

    std::vector<TensorType> u_; // disentagler
    std::vector<TensorType> w_; // isometries
    std::vector<TensorType> rho_; // density matrices
    std::vector<PairSizeType> vec_;
}; //class

} //namespace

#endif // MERALAYER_H
