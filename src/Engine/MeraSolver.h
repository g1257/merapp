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
	    : params_(params), srep_(""), meraLayer_(params.tauMax,0)
	{
		for (SizeType i = 0; i < params.tauMax; ++i) {
			SizeType sites = calcSitesForLayer(i);
			meraLayer_[i] = new MeraLayerType(params,i,sites);
			std::cout<<(*meraLayer_[i]);
		}
	}

	~MeraSolver()
	{
		for (SizeType i = 0; i < params_.tauMax; ++i) {
			delete meraLayer_[i];
			meraLayer_[i] = 0;
		}
	}

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

	SizeType calcSitesForLayer(SizeType tau) const
	{
		if (tau == 0) return params_.numberOfSites;
		return meraLayer_[tau - 1]->outputSites();
	}

private:

	MeraSolver(const MeraSolver&);

	MeraSolver& operator=(const MeraSolver&);

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

	const ParametersForSolver& params_;
	PsimagLite::String srep_;
    VectorMeraLayerType meraLayer_;
}; //class

} //namespace

#endif
