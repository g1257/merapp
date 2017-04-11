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
#include "TensorEvalSlow.h"
#include "TensorEvalNew.h"
#include "TensorOptimizer.h"
#include "InputCheck.h"
#include "ModelSelector.h"
#include "ModelBase.h"

namespace Mera {

template<typename ComplexOrRealType>
class MeraSolver {

	typedef PsimagLite::InputNg<InputCheck> InputNgType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef TensorOptimizer<ComplexOrRealType,InputNgType::Readable> TensorOptimizerType;
	typedef typename PsimagLite::Vector<TensorOptimizerType*>::Type VectorTensorOptimizerType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename TensorOptimizerType::MapPairStringSizeType MapPairStringSizeType;
	typedef typename TensorOptimizerType::ParametersForSolverType ParametersForSolverType;
	typedef typename TensorOptimizerType::SrepStatementType SrepStatementType;
	typedef typename TensorOptimizerType::VectorSrepStatementType VectorSrepStatementType;
	typedef typename TensorOptimizerType::TensorEvalBaseType TensorEvalBaseType;
	typedef typename TensorEvalBaseType::PairStringSizeType PairStringSizeType;
	typedef typename TensorEvalBaseType::VectorPairStringSizeType VectorPairStringSizeType;
	typedef typename TensorEvalBaseType::TensorType TensorType;
	typedef typename TensorEvalBaseType::VectorTensorType VectorTensorType;
	typedef ModelBase<ComplexOrRealType> ModelBaseType;
	typedef ModelSelector<ModelBaseType> ModelType;
	typedef typename TensorOptimizerType::SymmetryLocalType SymmetryLocalType;
	typedef typename TensorOptimizerType::ParametersForMeraType ParametersForMeraType;

	static const int EVAL_BREAKUP = TensorOptimizerType::EVAL_BREAKUP;

public:

	MeraSolver(PsimagLite::String filename)
	    : paramsForMera_(filename),
	      symmLocal_(filename),
	      isMeraPeriodic_(false),
	      iterMera_(1),
	      iterTensor_(1),
	      indexOfRootTensor_(0),
	      model_(paramsForMera_.model, paramsForMera_.hamiltonianConnection),
	      paramsForLanczos_(0),
	      noSymmLocal_(false)
	{
		InputCheck inputCheck;
		InputNgType::Writeable ioWriteable(filename,inputCheck);
		InputNgType::Readable io(ioWriteable);

		paramsForLanczos_ = new ParametersForSolverType(io,"Mera");

		int x = 0;
		io.readline(x, "IsMeraPeriodic=");
		isMeraPeriodic_ = (x > 0);

		try {
			io.readline(iterMera_,"IterMera=");
		} catch (std::exception&) {}

		try {
			io.readline(iterTensor_,"IterTensor=");
		} catch (std::exception&) {}


		try {
			io.readline(x,"NoSymmetryLocal=");
		} catch (std::exception&) {}

		noSymmLocal_ = (x == 0) ? false : true;

		SymmetryLocalType* symmLocal = (noSymmLocal_) ? 0 : &symmLocal_;

		PsimagLite::String dstr("");
		io.readline(dstr,"DimensionSrep=");
		TensorSrep tdstr(dstr);
		findTensors(tdstr);

		initTensorNameIds();

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
				SizeType ignoreTerm = terms + 1;
				io.readline(ignoreTerm,"IgnoreTerm=");
				energyTerms_.resize(terms,0);

				PsimagLite::String findStr = "Environ=";
				for (SizeType i = 0; i < terms; ++i) {
					PsimagLite::String srep;
					io.readline(srep,findStr);
					if (i == ignoreTerm) continue;
					size_t index = srep.find("equal");
					if (index != PsimagLite::String::npos)
						srep.replace(index,5,"=");
					energyTerms_[i] = new SrepStatementType(srep);
					allTensorsDefinedOrDie(energyTerms_[i]->rhs());
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
			                                                   *paramsForLanczos_,
			                                                   paramsForMera_,
			                                                   *symmLocal));

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

		for (SizeType i = 0; i < energyTerms_.size(); ++i) {
			delete energyTerms_[i];
			energyTerms_[i] = 0;
		}

		delete paramsForLanczos_;
	}

	void optimize()
	{
		RealType eprev = 1e6;
		for (SizeType i = 0; i < iterMera_; ++i)
			optimizeAllTensors(i, eprev);
	}

private:

	void optimizeAllTensors(SizeType iter, RealType& eprev)
	{
		static bool seenRoot = false;
		bool optimizeOnlyFirstOfLayer = isMeraPeriodic_;

		if (paramsForMera_.options.find("OptimizeAllLayers") != PsimagLite::String::npos)
			optimizeOnlyFirstOfLayer = false;

		SizeType ntensors = tensorOptimizer_.size();

		for (SizeType i = 0; i < ntensors; ++i) {

			PsimagLite::String name = tensorOptimizer_[i]->nameId().first;
			SizeType id = tensorOptimizer_[i]->nameId().second;

			SizeType firstOfLayer = tensorOptimizer_[i]->firstOfLayer();
			if (optimizeOnlyFirstOfLayer && firstOfLayer != id && name != "r") {
				tensorOptimizer_[i]->copyFirstOfLayer(name, firstOfLayer);
				continue;
			}

			if (!seenRoot) {
				if (name != "r")
					continue;
				else
					seenRoot = true;
			}

			tensorOptimizer_[i]->optimize(iterTensor_,
			                              iter,
			                              paramsForMera_.evaluator);

			RealType e = energy();
			if (e > eprev) {
				std::cerr<<"MeraSolver: found larger energy ";
				std::cerr<<e<<" restoring previous...\n";
				tensorOptimizer_[i]->restoreTensor();
				e = energy();
			}

			eprev = e;
			PsimagLite::String str("energy after optimizing ");
			str += name + ttos(id) + "= ";
			std::cout<<str<<e<<" [ Remember shift=";
			std::cout<<model_().energyShift()<<" ]\n";
		}
	}

	class ParallelEnergyHelper {

	public:

		ParallelEnergyHelper(bool noSymmLocal,
		                     SymmetryLocalType& symmLocal,
		                     VectorSrepStatementType& energyTerms,
		                     const VectorPairStringSizeType& tensorNameAndIds,
		                     MapPairStringSizeType& nameIdsTensor,
		                     VectorTensorType& tensors,
		                     const ParametersForMeraType& paramsForMera)
		    : noSymmLocal_(noSymmLocal),
		      symmLocal_(symmLocal),
		      energyTerms_(energyTerms),
		      tensorNameIds_(tensorNameAndIds),
		      nameIdsTensor_(nameIdsTensor),
		      tensors_(tensors),
		      paramsForMera_(paramsForMera),
		      e_(PsimagLite::Concurrency::npthreads, 0.0)
		{}

		void doTask(SizeType taskNumber, SizeType threadNum)
		{
			assert(threadNum < e_.size());
			e_[threadNum] += energy(taskNumber);
		}

		SizeType tasks() const { return energyTerms_.size(); }

		void sync()
		{
			if (e_.size() == 0) return;
			for (SizeType i = 1; i < e_.size(); ++i)
				e_[0] += e_[i];
		}

		RealType energy() const
		{
			if (e_.size() == 0)
				throw PsimagLite::RuntimeError("ParallelEnergyHelper failed\n");
			return e_[0];
		}

	private:

		RealType energy(SizeType ind)
		{
			SymmetryLocalType* symmLocal = (noSymmLocal_) ? 0 : &symmLocal_;
			assert(ind < energyTerms_.size());
			SrepStatementType* ptr = energyTerms_[ind];
			if (!ptr) return 0.0;
			TensorEvalBaseType* tensorEval =
			        TensorOptimizerType::getTensorEvalPtr(paramsForMera_.evaluator,
			                                              *ptr,
			                                              tensors_,
			                                              tensorNameIds_,
			                                              nameIdsTensor_,
			                                              *symmLocal);

			typename TensorEvalBaseType::HandleType handle = tensorEval->operator()();
			while (!handle.done());
			delete tensorEval;
			tensorEval = 0;
			VectorSizeType args(1,0);
			return tensors_[nameIdsTensor_[PairStringSizeType("e",ind)]]->operator()(args);
		}

		bool noSymmLocal_;
		SymmetryLocalType& symmLocal_;
		VectorSrepStatementType& energyTerms_;
		const VectorPairStringSizeType& tensorNameIds_;
		MapPairStringSizeType& nameIdsTensor_;
		VectorTensorType& tensors_;
		const ParametersForMeraType& paramsForMera_;
		VectorRealType e_;
	}; // class ParallelEnergyHelper

	RealType energy()
	{
		typedef PsimagLite::Parallelizer<ParallelEnergyHelper> ParallelizerType;
		ParallelizerType threadedEnergies(PsimagLite::Concurrency::npthreads,
		                                  PsimagLite::MPI::COMM_WORLD);

		ParallelEnergyHelper parallelEnergyHelper(noSymmLocal_,
		                                          symmLocal_,
		                                          energyTerms_,
		                                          tensorNameIds_,
	                                              nameIdsTensor_,
		                                          tensors_,
		                                          paramsForMera_);

		threadedEnergies.loopCreate(parallelEnergyHelper);
		parallelEnergyHelper.sync();
		return parallelEnergyHelper.energy();
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
				tensors_[ind]->setToMatrix(model_().twoSiteHam(id));
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

	void allTensorsDefinedOrDie(const TensorSrep& srep)
	{
		SizeType ntensors = srep.size();
		for (SizeType i = 0; i < ntensors; ++i) {
			PsimagLite::String name = srep(i).name();
			SizeType id = srep(i).id();
			PairStringSizeType p(name,id);
			SizeType x = nameIdsTensor_[p];
			if (tensorNameIds_[x] == p) continue;
			PsimagLite::String msg("Not found tensor "+ name + ttos(id) + "\n");
			throw PsimagLite::RuntimeError("allTensorsDefinedOrDie: " + msg);
		}
	}

	MeraSolver(const MeraSolver&);

	MeraSolver& operator=(const MeraSolver&);

	const ParametersForMeraType paramsForMera_;
	SymmetryLocalType symmLocal_;
	bool isMeraPeriodic_;
	SizeType iterMera_;
	SizeType iterTensor_;
	SizeType indexOfRootTensor_;
	VectorPairStringSizeType tensorNameIds_;
	MapPairStringSizeType nameIdsTensor_;
	VectorTensorType tensors_;
	ModelType model_;
	VectorTensorOptimizerType tensorOptimizer_;
	ParametersForSolverType* paramsForLanczos_;
	bool noSymmLocal_;
	VectorSrepStatementType energyTerms_;
}; // class MeraSolver
} // namespace Mera
#endif // MERASOLVER_H
