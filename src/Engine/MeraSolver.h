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
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	MeraSolver(const ParametersForSolver& params)
	    : params_(params), srep_(""), meraLayer_(params.tauMax,0)
	{
		for (SizeType i = 0; i < params.tauMax; ++i) {
			SizeType sites = calcSitesForLayer(i);
			meraLayer_[i] = new MeraLayerType(params,i,sites,(i > 0) ? meraLayer_[i-1] : 0);
			std::cout<<(*meraLayer_[i]);
			srep_ += meraLayer_[i]->sRep();
		}

		cleanUnpaired(srep_);
	}

	~MeraSolver()
	{
		for (SizeType i = 0; i < params_.tauMax; ++i) {
			delete meraLayer_[i];
			meraLayer_[i] = 0;
		}
	}

	const PsimagLite::String& sRep() const { return srep_; }

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

	friend std::ostream& operator<<(std::ostream& os, const MeraSolver& ms)
	{
		os<<"fullSrep="<<ms.srep_<<"\n";
		return os;
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

	void cleanUnpaired(PsimagLite::String& srep)
	{
		SizeType max = findMaxSummed(srep);
		VectorSizeType counter(max,0);
		getSpairCount(counter,srep);

		for (SizeType i = 0; i < max; ++i) {
			if (counter[i] == 2) continue;
			std::cerr<<"Upaired s"<<i<<" with "<<counter[i]<<"\n";
		}

		PsimagLite::String srep2;
		SizeType i = 0;
		while (i < srep.length()) {
			if (srep[i] == 's') {
				SizeType j = i + 1;
				PsimagLite::String ds("");
				while (j < srep.length()) {
					if (!isdigit(srep[j])) break;
					ds += srep[j];
					++j;
				}

				SizeType d = atoi(ds.c_str());
				assert(d < counter.size());
				if (counter[d] == 2) {
					srep2 += "s";
					++i;
					continue;
				}

				srep2 += "d";
				i = j - 1;
			}

			srep2 += srep[i];
			++i;
		}

		srep = srep2;
	}

	void getSpairCount(VectorSizeType& counter,PsimagLite::String srep) const
	{
		SizeType i = 0;
		while (i < srep.length()) {
			if (srep[i] == 's') {
				SizeType j = i + 1;
				PsimagLite::String ds("");
				while (j < srep.length()) {
					if (!isdigit(srep[j])) break;
					ds += srep[j];
					++j;
				}

				i = j;
				SizeType d = atoi(ds.c_str());
				assert(d < counter.size());
				counter[d]++;
				continue;
			}

			++i;
		}
	}

	SizeType findMaxSummed(PsimagLite::String& srep)
	{
		SizeType i = 0;
		SizeType max = 0;
		while (i < srep.length()) {
			if (srep[i] == 's') {
				SizeType j = i + 1;
				PsimagLite::String ds("");
				while (j < srep.length()) {
					if (!isdigit(srep[j])) break;
					ds += srep[j];
					++j;
				}

				i = j;
				SizeType d = atoi(ds.c_str());
				if (d > max) max = d;
				continue;
			}

			++i;
		}

		return max + 1;
	}

	const ParametersForSolver& params_;
	PsimagLite::String srep_;
    VectorMeraLayerType meraLayer_;
}; //class

} //namespace

#endif
