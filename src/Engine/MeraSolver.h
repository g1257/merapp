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
#include "DimensionSrep.h"

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
	typedef TensorEvalSlow<ComplexOrRealType> TensorEvalSlowType;
	typedef TensorEvalNew<ComplexOrRealType> TensorEvalNewType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	static const int EVAL_BREAKUP = TensorOptimizerType::EVAL_BREAKUP;

public:

	MeraSolver(PsimagLite::String filename)
	    : paramsForMera_(filename),
	      symmLocal_(0),
	      meraStr_(""),
	      isMeraPeriodic_(false),
	      iterMera_(1),
	      iterTensor_(1),
	      indexOfRootTensor_(0),
	      nameToIndexLut_(0),
	      model_(paramsForMera_.model, paramsForMera_.hamiltonianConnection),
	      paramsForLanczos_(0)
	{
		TensorType::init();
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


		io.readline(meraStr_,"MERA=");
		io.readline(dsrepEnvirons_, "DsrepEnvirons=");

		PsimagLite::String hString = "D" + ttos(model_().qOne().size());
		PsimagLite::String args = "(" + hString + "," + hString + "|";
		args += hString + "," + hString + ")";
		for (SizeType i = 0; i < paramsForMera_.hamiltonianConnection.size(); ++i) {
			if (paramsForMera_.hamiltonianConnection[i] == 0.0) continue;
			meraStr_ += "h" + ttos(i) + args;
		}

		try {
			io.readline(x,"NoSymmetryLocal=");
		} catch (std::exception&) {}

		bool noSymmLocal = (x == 0) ? false : true;

		updateTensorSizes(noSymmLocal);

		nameToIndexLut_ = new NameToIndexLut<TensorType>(tensors_);
		bool rootTensorFound = false;
		while (true) {
			PsimagLite::String str("");
			try {
				io.readline(str,"TensorId=");
			} catch (std::exception&) {
				break;
			}

			PsimagLite::Vector<PsimagLite::String>::Type tokens;
			PsimagLite::split(tokens, str, ",");
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
			                                                   tensors_,
			                                                   *nameToIndexLut_,
			                                                   *paramsForLanczos_,
			                                                   paramsForMera_,
			                                                   symmLocal_));

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

		TensorType::finalize();
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

		ParallelEnergyHelper(SymmetryLocalType* symmLocal,
		                     VectorSrepStatementType& energyTerms,
		                     VectorTensorType& tensors,
		                     NameToIndexLut<TensorType>& nameToIndexLut,
		                     const ParametersForMeraType& paramsForMera)
		    : symmLocal_(symmLocal),
		      energyTerms_(energyTerms),
		      tensors_(tensors),
		      nameToIndexLut_(nameToIndexLut),
		      paramsForMera_(paramsForMera),
		      e_(PsimagLite::Concurrency::codeSectionParams.npthreads, 0.0)
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
			assert(ind < energyTerms_.size());
			SrepStatementType* ptr = energyTerms_[ind];
			if (!ptr) return 0.0;
			TensorEvalBaseType* tensorEval = TensorOptimizerType::ParallelEnvironHelperType::
			        getTensorEvalPtr(paramsForMera_.evaluator,
			                         *ptr,
			                         tensors_,
			                         nameToIndexLut_,
			                         symmLocal_);

			typename TensorEvalBaseType::HandleType handle = tensorEval->operator()();
			while (!handle.done());

			VectorSizeType args(1,0);
			const SizeType index = tensorEval->nameToIndexLut("e" + ttos(ind));
			delete tensorEval;
			tensorEval = 0;
			assert(index < tensors_.size());
			return tensors_[index]->operator()(args);
		}

		SymmetryLocalType* symmLocal_;
		VectorSrepStatementType& energyTerms_;
		VectorTensorType& tensors_;
		NameToIndexLut<TensorType>& nameToIndexLut_;
		const ParametersForMeraType& paramsForMera_;
		VectorRealType e_;
	}; // class ParallelEnergyHelper

	RealType energy()
	{
		typedef PsimagLite::Parallelizer<ParallelEnergyHelper> ParallelizerType;
		ParallelizerType threadedEnergies(PsimagLite::Concurrency::codeSectionParams);

		ParallelEnergyHelper parallelEnergyHelper(symmLocal_,
		                                          energyTerms_,
		                                          tensors_,
		                                          *nameToIndexLut_,
		                                          paramsForMera_);

		threadedEnergies.loopCreate(parallelEnergyHelper);
		parallelEnergyHelper.sync();
		return parallelEnergyHelper.energy();
	}

	void initTensorNameIds(VectorStringType& tensorNameIds)
	{
		PsimagLite::Sort<VectorStringType> sort;
		VectorSizeType perm(tensorNameIds.size(),0);
		sort.sort(tensorNameIds,perm);
		SizeType end = (std::unique(tensorNameIds.begin(),
		                            tensorNameIds.end()) - tensorNameIds.begin());
		tensorNameIds.resize(end);
	}

	void initTensors(const VectorStringType& tensorNameIds,
	                 const TensorSrep& td)
	{
		tensors_.resize(tensorNameIds.size());
		SizeType ntensors = tensors_.size();

		if (td.size() != ntensors) {
			PsimagLite::String str("TensorOptimizer dimension string " + ttos(td.size()));
			str += ", was expecting " + ttos(ntensors) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		for (SizeType i = 0; i < ntensors; ++i) {
			TensorSrep::PairStringSizeType mypair = TensorSrep::
			        splitIntoNameAndId(td(i).fullName());
			typename VectorStringType::const_iterator x = std::find(tensorNameIds.begin(),
			                                                        tensorNameIds.end(),
			                                                        td(i).fullName());
			if (x == tensorNameIds.end()) {
				std::cerr<<"WARNING: Unused tensor fullname= "<< td(i).fullName()<<"\n";
				continue;
			}

			SizeType ind = x - tensorNameIds.begin();
			assert(ind < tensors_.size());

			SizeType legs = td(i).legs();
			VectorSizeType dimensions(legs,0);
			for (SizeType j = 0; j < legs; ++j) {
				SizeType legTag = td(i).legTag(j);
				dimensions[j] = legTag;
			}

			SizeType ins = td(i).ins();
			assert(ind < tensors_.size());
			tensors_[ind] = new TensorType(td(i).fullName(), dimensions, ins);
			if (td(i).fullName()[0] == 'h') {
				tensors_[ind]->setToMatrix(model_().twoSiteHam(mypair.second));
			} else {
				tensors_[ind]->setToIdentity(1.0);
			}
		}
	}

	void findTensors(VectorStringType& tensorNameIds,
	                 const TensorSrep& t)
	{
		SizeType ntensors = t.size();
		for (SizeType i = 0; i < ntensors; ++i)
			tensorNameIds.push_back(t(i).fullName());
	}

	void allTensorsDefinedOrDie(const TensorSrep& srep)
	{
		SizeType ntensors = srep.size();
		for (SizeType i = 0; i < ntensors; ++i)
			nameToIndexLut_->operator()(srep(i).fullName());
	}

	void updateTensorSizes(bool noSymmLocal)
	{
		TensorSrep tsrep(meraStr_);
		SizeType maxLegs = 2.0*paramsForMera_.hamiltonianConnection.size();
		SymmetryLocal* symmLocal = new SymmetryLocal(tsrep.size(), model_().qOne(), maxLegs);
		DimensionSrep<SymmetryLocal> dimSrep(meraStr_, *symmLocal, paramsForMera_.m);
		PsimagLite::String dsrep = dimSrep() + dsrepEnvirons_;

		delete symmLocal_;
		symmLocal_ = 0;
		if (!noSymmLocal)
			symmLocal_ = symmLocal;

		TensorSrep tdstr(dsrep);
		VectorStringType tensorNameIds;
		findTensors(tensorNameIds,  tdstr);
		initTensorNameIds(tensorNameIds);

		initTensors(tensorNameIds, tdstr);
	}

	MeraSolver(const MeraSolver&);

	MeraSolver& operator=(const MeraSolver&);

	const ParametersForMeraType paramsForMera_;
	SymmetryLocalType* symmLocal_;
	PsimagLite::String meraStr_;
	PsimagLite::String dsrepEnvirons_;
	bool isMeraPeriodic_;
	SizeType iterMera_;
	SizeType iterTensor_;
	SizeType indexOfRootTensor_;
	VectorTensorType tensors_;
	NameToIndexLut<TensorType>* nameToIndexLut_;
	ModelType model_;
	VectorTensorOptimizerType tensorOptimizer_;
	ParametersForSolverType* paramsForLanczos_;
	VectorSrepStatementType energyTerms_;
}; // class MeraSolver
} // namespace Mera
#endif // MERASOLVER_H
