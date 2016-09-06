#ifndef MERASOLVER_H
#define MERASOLVER_H
#include <iostream>
#include "MeraLayer.h"
#include "ParametersForSolver.h"

namespace Mera {

class MeraSolver {

	typedef ParametersForSolver ParametersForSolverType;
	typedef MeraLayer<ParametersForSolverType> MeraLayerType;
	typedef typename PsimagLite::Vector<MeraLayerType*>::Type VectorMeraLayerType;
    typedef typename MeraLayerType::PairSizeType PairSizeType;

public:

	MeraSolver(const ParametersForSolver& params)
	{}

    void computeGroundState()
    {
        // Compute All Density Matricies - top to bottom using A or D operators
        SizeType qiter = 1;
        SizeType layers = 1;

        for (SizeType iter=0; iter<qiter; ++iter) {
            for (SizeType layer=0; layer<layers; ++layer) {
                optimize(iter,layer);
            }
        }


    }

private:

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
