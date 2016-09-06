#ifndef MERASTEP_H
#define MERASTEP_H
#include "Vector.h"
#include "MeraLayer.h"
#include <cassert>

namespace Mera {

template<typename ParametersForSolverType>
class MeraStep {

    typedef MeraLayer MeraLayerType;
	typedef typename PsimagLite::Vector<MeraLayerType*>::Type VectorMeraLayerType;
    typedef MeraLayerType::PairSizeType PairSizeType;

public:

	MeraStep(const ParametersForSolverType& params)
	    : meraLayer_(params.tauMax)
	{

	}

    void optimize(SizeType iter, SizeType layer)
    {
        SizeType qlayer = 1; // number of loops for convergece of u or w
        SizeType n = meraLayer_[layer]->size();
        for (SizeType i=0; i<qlayer; ++i) {
            // optimize one (w or u) of a layer tau
            for (SizeType f=0; f<n; ++f) {
                // factors = number of u and w in a layer
                optimizeThisTensor(meraLayer_[layer]->tensorToOptimize(f));
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

    VectorMeraLayerType meraLayer_;

}; //class

} //namespace

#endif
