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
#include "InputNg.h"
#include "TensorSrep.h"
#include "TensorEval.h"
#include "TensorOptimizer.h"
#include "InputCheck.h"

namespace Mera {

template<typename ComplexOrRealType>
class MeraSolver {

	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef TensorEval<ComplexOrRealType> TensorEvalType;
	typedef typename TensorEvalType::TensorType TensorType;
	typedef TensorOptimizer<ComplexOrRealType,InputNgType::Readable> TensorOptimizerType;
	typedef typename PsimagLite::Vector<TensorOptimizerType*>::Type VectorTensorOptimizerType;
	typedef typename TensorEvalType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalType::VectorTensorType VectorTensorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename TensorOptimizerType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename TensorOptimizerType::ParametersForSolverType ParametersForSolverType;
	typedef typename TensorOptimizerType::SrepEquationType SrepEquationType;
	typedef typename TensorOptimizerType::VectorSrepEquationType VectorSrepEquationType;

	static const int EVAL_BREAKUP = TensorOptimizerType::EVAL_BREAKUP;

public:

	MeraSolver(PsimagLite::String filename)
	    : iterMera_(1),
	      iterTensor_(1),
	      indexOfRootTensor_(0),
	      twoSiteHam_(4,4),
	      paramsForLanczos_(0)
	{
		InputCheck inputCheck;
		InputNgType::Writeable ioWriteable(filename,inputCheck);
		InputNgType::Readable io(ioWriteable);

		paramsForLanczos_ = new ParametersForSolverType(io,"Mera");

		try {
			io.readline(iterMera_,"IterMera=");
		} catch (std::exception&) {}

		try {
			io.readline(iterTensor_,"IterTensor=");
		} catch (std::exception&) {}

		PsimagLite::String dstr("");
		io.readline(dstr,"DimensionSrep=");
		TensorSrep tdstr(dstr);
		findTensors(tdstr);

		initTensorNameIds();

		bool makeHamTheIdentity = false;
		setTwoSiteHam(makeHamTheIdentity);

		initTensors(tdstr);

		bool rootTensorFound = false;
		while (true) {
			PsimagLite::String str("");
			try {
				io.readline(str,"TensorId=");
			} catch (std::exception&) {
				break;
			}

			PsimagLite::Vector<PsimagLite::String>::Type tokens;
			PsimagLite::tokenizer(str,tokens,",");
			if (tokens.size() != 2) {
				PsimagLite::String str("MeraSolver: Error reading TensorId=");
				str += "from file " + filename + "\n";
				throw PsimagLite::RuntimeError(str);
			}

			PsimagLite::String name = tokens[0];
			if (name == "E") {
				SizeType terms = 0;
				io.readline(terms,"Terms=");

				energyTerms_.resize(terms,0);

				PsimagLite::String findStr = "Environ=";
				for (SizeType i = 0; i < terms; ++i) {
					PsimagLite::String srep;
					io.readline(srep,findStr);
					size_t index = srep.find("equal");
					if (index != PsimagLite::String::npos)
						srep.replace(index,5,"=");
					energyTerms_[i] = new SrepEquationType(srep);
				}

				continue;
			}

			if (name != "u" && name != "w" && name != "r")
				continue;

			SizeType id = atoi(tokens[1].c_str());
			tensorOptimizer_.push_back(new TensorOptimizerType(io,
			                                                   name,
			                                                   id,
			                                                   tensorNameIds_,
			                                                   nameIdsTensor_,
			                                                   tensors_,
			                                                   *paramsForLanczos_));

			if (name == "r") {
				if (rootTensorFound) {
					PsimagLite::String msg("FATAL: File " + filename);
					throw PsimagLite::RuntimeError(msg + " more than one root found\n");
				}

				assert(tensorOptimizer_.size() > 0);
				indexOfRootTensor_ = tensorOptimizer_.size() - 1;
				rootTensorFound = true;
			}
		}

		if (!rootTensorFound) {
			PsimagLite::String msg("FATAL: File " + filename);
			throw PsimagLite::RuntimeError(msg + " root tensor not found\n");
		}

		if (energyTerms_.size() == 0) {
			PsimagLite::String msg("FATAL: File " + filename);
			throw PsimagLite::RuntimeError(msg + " energyTerms not found\n");
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

		delete paramsForLanczos_;
	}

	void optimize()
	{
		for (SizeType i = 0; i < iterMera_; ++i)
			optimizeAllTensors(i);
	}

private:

	void optimizeAllTensors(SizeType iter)
	{
		SizeType ntensors = tensorOptimizer_.size();
//		SizeType prevLayer = 0;
		for (SizeType i = 0; i < ntensors; ++i) {
//			SizeType thisLayer = tensorOptimizer_[i]->layer();
//			if (prevLayer != thisLayer || i == 0) {
//				tensorOptimizer_[indexOfRootTensor_]->optimize(iterTensor_,
//				                                               iter);
//				prevLayer = thisLayer;
//			}

//			if (i != indexOfRootTensor_)
				tensorOptimizer_[i]->optimize(iterTensor_,iter);

			PsimagLite::String str("energy after optimizing ");
			str += tensorOptimizer_[i]->nameId().first;
			str += ttos(tensorOptimizer_[i]->nameId().second) + "= ";
			std::cout<<str<<energy()<<"\n";
		}
	}

	RealType energy()
	{
		RealType e = 0;
		for (SizeType i = 0; i < energyTerms_.size(); ++i)
			e += energy(i);

		return e;
	}

	RealType energy(SizeType ind)
	{
		TensorEvalType tensorEval(*(energyTerms_[ind]),
		                          tensors_,
		                          tensorNameIds_,
		                          nameIdsTensor_,
		                          EVAL_BREAKUP);
		typename TensorEvalType::HandleType handle = tensorEval(EVAL_BREAKUP);
		while (!handle.done());
		VectorSizeType args(1,0);
		return tensors_[nameIdsTensor_[PairStringSizeType("e",ind)]]->operator()(args);
	}

	// FIXME: pick up model dependency here
	void setTwoSiteHam(bool testWithIdentity)
	{
		SizeType n = twoSiteHam_.n_row();
		assert(n = twoSiteHam_.n_col());
		for (SizeType i = 0; i < n; ++i)
			twoSiteHam_(i,i) = 1.0;
		if (testWithIdentity) return;

		// Sz Sz
		for (SizeType i = 0; i < n; ++i)
			twoSiteHam_(i,i) = (i == 0 || i ==3) ? 0.25 : -0.25;

		// S+S- S-S+
		twoSiteHam_(1,2) = twoSiteHam_(2,1) = 0.5;

		normalizeHam();
	}

	void normalizeHam()
	{
		MatrixType m = twoSiteHam_;
		SizeType n = m.n_row();
		VectorRealType eigs(n,0.0);
		diag(m,eigs,'N');
		assert(n - 1 < eigs.size());
		RealType diagCorrection = eigs[n-1];
		std::cout<<"MeraSolver: DiagonalCorrection= "<<diagCorrection<<"\n";
		for (SizeType i = 0; i < n; ++i)
			twoSiteHam_(i,i) -= diagCorrection;
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
		for (SizeType i = 0; i < tensorNameIds_.size(); ++i)
			nameIdsTensor_[tensorNameIds_[i]] = i;
	}

	void initTensors(const TensorSrep& td)
	{
		tensors_.resize(tensorNameIds_.size());
		SizeType ntensors = tensors_.size();

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

			SizeType legs = td(i).legs();
			VectorSizeType dimensions(legs,0);
			for (SizeType j = 0; j < legs; ++j) {
				SizeType legTag = td(i).legTag(j);
				dimensions[j] = legTag;
			}

			SizeType ins = td(i).ins();
			assert(ind < tensors_.size());
			tensors_[ind] = new TensorType(dimensions,ins);
			if (name == "h") {
				tensors_[ind]->setToMatrix(twoSiteHam_);
			} else if (name == "r") {
				assert(0 < dimensions.size());
				tensors_[ind]->setToIdentity(1.0);
			} else if (name == "u") {
				tensors_[ind]->setToIdentity(1.0);
			} else {
				tensors_[ind]->setToIdentity(1.0);
			}
		}
	}

	void findTensors(const TensorSrep& t)
	{
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

	SizeType iterMera_;
	SizeType iterTensor_;
	SizeType indexOfRootTensor_;
	VectorPairStringSizeType tensorNameIds_;
	MapPairStringSizeType nameIdsTensor_;
	VectorTensorType tensors_;
	MatrixType twoSiteHam_;
	VectorTensorOptimizerType tensorOptimizer_;
	ParametersForSolverType* paramsForLanczos_;
	VectorSrepEquationType energyTerms_;
}; // class MeraSolver
} // namespace Mera
#endif // MERASOLVER_H
