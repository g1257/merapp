#ifndef MERASTEP_H
#define MERASTEP_H
#include "Vector.h"
#include "MeraLayer.h"

namespace Mera {

class MeraStep {

    typedef MeraLayer MeraLayerType;
    typedef std::pair<SizeType,SizeType> PairSizeType;

    enum TensorTypeEnum {TENSOR_TYPE_W,TENSOR_TYPE_U};
public:

    void optimize(SizeType iter, SizeType layer)
    {
        SizeType n = meraLayer_.size();
        SizeType nOverTwo = n/2;
        std::vector<PairSizeType> vec(n);

        for (SizeType i=0; i<nOverTwo; ++i) {
            vec[2*i] = PairSizeType(i,TENSOR_TYPE_W);
            vec[2*i+1] = PairSizeType(i,TENSOR_TYPE_U);
        }

        SizeType qlayer = 1; // number of loops for convergece of u or w
        for (SizeType i=0; i<qlayer; ++i) {
            // optimize one (w or u) of a layer tau
            for (SizeType f=0; f<n; ++f) {
                // factors = number of u and w in a layer
                optimizeThisTensor(vec[f]);
            }
        }
    }

private:

    void optimizeThisTensor(PairSizeType pair)
    {
        // find Y (environment) for this tensor
        // Y = USV^+
        // w = -VU^+
    }

    MeraLayerType meraLayer_;

}; //class

} //namespace

#endif
