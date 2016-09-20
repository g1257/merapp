/*
Copyright (c) 2016, UT-Battelle, LLC

MERA++, Version 0.

This file is part of MERA++.
MERA++ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MERA++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MERA++. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MERASOLVER_H
#define MERASOLVER_H
#include "IoSimple.h"
#include "TensorSrep.h"
#include "TensorEval.h"
#include "TensorOptimizer.h"

namespace Mera {

template<typename ComplexOrRealType>
class MeraSolver {

	typedef PsimagLite::IoSimple::In IoInType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef TensorEval<ComplexOrRealType> TensorEvalType;
	typedef typename TensorEvalType::TensorType TensorType;
	typedef TensorOptimizer<ComplexOrRealType> TensorOptimizerType;
	typedef typename PsimagLite::Vector<TensorOptimizerType*>::Type VectorTensorOptimizerType;
	typedef typename TensorEvalType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalType::VectorTensorType VectorTensorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

public:

	MeraSolver(PsimagLite::String file)
	    : tauMax_(0), iterMera_(2), iterTensor_(5), twoSiteHam_(4,4)
	{
		IoInType io(file);
		io.readline(tauMax_,"TauMax=");

		try {
			io.readline(iterMera_,"IterMera=");
		} catch (std::exception&) {
			io.rewind();
		}

		try {
			io.readline(iterTensor_,"IterTensor=");
		} catch (std::exception&) {
			io.rewind();
		}

		PsimagLite::String srepMera;
		io.readline(srepMera,"MERA=");
		srepMera += "h0(2,2,2,2)";
		findTensors(srepMera);

		bool makeHamTheIdentity = false;
		setTwoSiteHam(makeHamTheIdentity);

		initTensorNameIds();

		PsimagLite::String dstr("");
		io.rewind();
		io.readline(dstr,"DIMENSION_SREP=");
		initTensors(dstr);

		while (!io.eof()) {
			PsimagLite::String str("");
			try {
				io.readline(str,"Y_FOR_TENSOR=");
			} catch (std::exception&) {
				break;
			}

			PsimagLite::Vector<PsimagLite::String>::Type tokens;
			PsimagLite::tokenizer(str,tokens,",");
			if (tokens.size() != 2) {
				PsimagLite::String str("MeraSolver: Error reading Y_FOR_TENSOR");
				str += "from file " + file + "\n";
				throw PsimagLite::RuntimeError(str);
			}

			PsimagLite::String name = tokens[0];
			SizeType id = atoi(tokens[1].c_str());
			tensorOptimizer_.push_back(new TensorOptimizerType(io,
			                                                   name,
			                                                   id,
			                                                   tensorNameIds_,
			                                                   tensors_));
		}

		std::cerr<<"MeraSolver::ctor() done\n";
	}

	~MeraSolver()
	{
		for (SizeType i = 0; i < tensorOptimizer_.size(); ++i) {
			delete tensorOptimizer_[i];
			tensorOptimizer_[i] = 0;
		}

		for (SizeType i = 0; i < tensors_.size(); ++i) {
			delete tensors_[i];
			tensors_[i] = 0;
		}
	}

	void optimize()
	{
		for (SizeType i = 0; i < iterMera_; ++i)
			optimizeAllTensors();
	}

private:

	void optimizeAllTensors()
	{
		SizeType ntensors = tensorOptimizer_.size();
		for (SizeType i = 0; i < ntensors; ++i)
			tensorOptimizer_[i]->optimize(iterTensor_);
	}

	// FIXME: pick up model dependency here
	void setTwoSiteHam(bool testWithIdentity)
	{
		for (SizeType i = 0; i < 4; ++i)
			twoSiteHam_(i,i) = 1.0;
		if (testWithIdentity) return;

		// Sz Sz
		for (SizeType i = 0; i < 4; ++i)
			twoSiteHam_(i,i) = (i == 0 || i ==3) ? 0.25 : -0.25;

		// S+S- S-S+
		twoSiteHam_(1,2) = twoSiteHam_(2,1) = 0.5;
	}

	void initTensorNameIds()
	{
		PsimagLite::Sort<VectorPairStringSizeType> sort;
		VectorSizeType perm(tensorNameIds_.size(),0);
		sort.sort(tensorNameIds_,perm);
		SizeType end = (std::unique(tensorNameIds_.begin(),
		                            tensorNameIds_.end()) -
		                tensorNameIds_.begin());
		tensorNameIds_.resize(end);
	}

	void initTensors(PsimagLite::String dstr)
	{
		tensors_.resize(tensorNameIds_.size());
		SizeType ntensors = tensors_.size();

		TensorSrep td(dstr);
		if (td.size() != ntensors) {
			PsimagLite::String str("TensorOptimizer dimension string " + ttos(td.size()));
			str += ", was expecting " + ttos(ntensors) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		for (SizeType i = 0; i < ntensors; ++i) {
			PsimagLite::String name = td(i).name();
			SizeType id = td(i).id();
			PairStringSizeType p(name,id);
			typename VectorPairStringSizeType::iterator x = std::find(tensorNameIds_.begin(),
			                                                          tensorNameIds_.end(),
			                                                          p);
			if (x == tensorNameIds_.end()) {
				std::cerr<<"WARNING: Unused tensor name= "<<name<<" id= "<<id<<"\n";
				continue;
			}

			SizeType ind = x - tensorNameIds_.begin();
			assert(ind < tensors_.size());

			SizeType ins = td(i).ins();
			SizeType outs = td(i).outs();
			VectorSizeType dimensions(ins + outs);
			for (SizeType j = 0; j < ins; ++j) {
				SizeType legTag = td(i).legTag(j,TensorStanza::INDEX_DIR_IN);
				dimensions[j] = legTag;
			}

			for (SizeType j = 0; j < outs; ++j) {
				SizeType legTag = td(i).legTag(j,TensorStanza::INDEX_DIR_OUT);
				dimensions[j+ins] = legTag;
			}

			assert(ind < tensors_.size());
			tensors_[ind] = new TensorType(dimensions,ins);
			if (name == "h") {
				tensors_[ind]->setToMatrix(twoSiteHam_);
			} else if (name == "r") {
				assert(0 < dimensions.size());
				tensors_[ind]->setToIdentityEx(1.0/sqrt(2.0),1);
			} else {
				tensors_[ind]->setToIdentity(1.0);
			}
		}
	}

	void findTensors(PsimagLite::String srep)
	{
		TensorSrep t(srep);
		SizeType ntensors = t.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			PsimagLite::String name = t(i).name();
			SizeType id = t(i).id();
			PairStringSizeType p(name,id);
			tensorNameIds_.push_back(p);
		}
	}

	MeraSolver(const MeraSolver&);

	MeraSolver& operator=(const MeraSolver&);

	SizeType tauMax_;
	SizeType iterMera_;
	SizeType iterTensor_;
	VectorPairStringSizeType tensorNameIds_;
	VectorTensorType tensors_;
	MatrixType twoSiteHam_;
	VectorTensorOptimizerType tensorOptimizer_;
}; // class MeraSolver
} // namespace Mera
#endif // MERASOLVER_H
