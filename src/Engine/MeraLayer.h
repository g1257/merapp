#ifndef MERALAYER_H
#define MERALAYER_H
#include <cassert>
#include "Vector.h"
#include "Matrix.h"
#include "TensorLeg.h"

namespace Mera {

template<typename ParametersForSolverType>
class MeraLayer {

	typedef int TensorType;
	typedef TensorLeg TensorLegType_;
	typedef PsimagLite::Matrix<TensorLegType_*> MatrixTensorLegType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ParametersForSolverType::MeraArchitectureEnum MeraArchitectureEnum;

	enum TensorTypeEnum {TENSOR_TYPE_W,TENSOR_TYPE_U};

public:

	typedef TensorLegType_ TensorLegType;
	typedef std::pair<SizeType,TensorTypeEnum> PairSizeTypeType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

	MeraLayer(const ParametersForSolverType& params,
	          SizeType tau,
	          SizeType sitesInLayer,
	          MeraLayer* prevLayer)
	    : params_(params),
	      tau_(tau),
	      sitesInLayer_(sitesInLayer),
	      prevLayer_(prevLayer),
	      outputSites_(0),
	      scount_(0),
	      srep_("")
	{
		setMeraArchitecture(params.meraArch);
		setUpdateOrder();
		createSrep(srep_);
	}

	~MeraLayer()
	{
		for (SizeType i = 0; i < mw_.n_row(); ++i) {
			for (SizeType j = 0; j < mw_.n_col(); ++j) {
				delete mw_(i,j);
				mw_(i,j) = 0;
			}
		}

		for (SizeType i = 0; i < mu_.n_row(); ++i) {
			for (SizeType j = 0; j < mu_.n_col(); ++j) {
				delete mu_(i,j);
				mu_(i,j) = 0;
			}
		}
	}

	SizeType size() const
	{
		return vecForUpdateOrder_.size();
	}

	SizeType outputSites() const { return outputSites_; }

	const PsimagLite::String sRep() const { return srep_; }

	const PairSizeTypeType& tensorToOptimize(SizeType i) const
	{
		assert(i < vecForUpdateOrder_.size());
		return vecForUpdateOrder_[i];
	}

	template<typename P>
	friend std::ostream& operator<<(std::ostream& os, const MeraLayer<P>& m)
	{
		os<<"tau="<<m.tau_<<"\n";
		os<<"sitesInThisLayer="<<m.sitesInLayer_<<"\n";
		os<<"vs="<<m.mw_.n_row()<<"\n";
		for (SizeType i = 0; i < m.mw_.n_row(); ++i) {
			os<<"isometry number "<<i<<" has sites ";
			for (SizeType j = 0; j < m.mw_.n_col(); ++j) {
				TensorLegType* ptr = m.mw_(i,j);
				if (!ptr) continue;
				PsimagLite::String val = (ptr->inOrOut == TensorLegType::OUT) ? " *" : " ";
				if (ptr->site >= m.sitesInLayer_)
					os<<val<<"o";
				else
					os<<val<<ptr->site;
			}

			os<<"\n";
		}

		os<<"us="<<m.mu_.n_row()<<"\n";
		for (SizeType i = 0; i < m.mu_.n_row(); ++i) {
			os<<"disentangler number "<<i<<" has sites ";
			for (SizeType j = 0; j < m.mu_.n_col(); ++j) {
				TensorLegType* ptr = m.mu_(i,j);
				if (!ptr) continue;
				PsimagLite::String val = (ptr->inOrOut == TensorLegType::OUT) ? " *" : " ";
				if (ptr->site >= m.sitesInLayer_)
					os<<val<<"o";
				else
					os<<val<<ptr->site;
			}

			os<<"\n";
		}

		os<<"srep="<<m.srep_<<"\n";
		os<<"-------------------\n";
		return os;
	}

private:

	MeraLayer(const MeraLayer&);

	MeraLayer& operator=(const MeraLayer&);

	void setUpdateOrder()
	{
		SizeType n = std::min(mu_.n_row(), mw_.n_row());
		vecForUpdateOrder_.resize(mu_.n_row() + mw_.n_row());

		for (SizeType i=0; i<n; ++i) {
			vecForUpdateOrder_[2*i] = PairSizeTypeType(i,TENSOR_TYPE_U);
			vecForUpdateOrder_[2*i+1] = PairSizeTypeType(i,TENSOR_TYPE_W);
		}

		int x = mu_.n_row() - mw_.n_row();
		if (x == 0) return;

		SizeType m = std::abs(x);
		TensorTypeEnum type = (x < 0) ? TENSOR_TYPE_W : TENSOR_TYPE_U;
		for (SizeType i=n; i<n+m; ++i)
			vecForUpdateOrder_[n+i] =  PairSizeTypeType(i,type);
	}

	void setMeraArchitecture(MeraArchitectureEnum what)
	{
		if (what == ParametersForSolverType::MERA_ARCH_1D_BINARY)
			setMeraArchitecture1dBinary();
		else if (what == ParametersForSolverType::MERA_ARCH_1D_TERNARY)
			setMeraArchitecture1dTernary();
		else
			throw PsimagLite::RuntimeError("MERA architecture not implemented\n");
		markBogusUoutputs();
	}

	void setMeraArchitecture1dBinary()
	{
		SizeType totalUs = (sitesInLayer_&1) ? (sitesInLayer_ - 1)/2 : sitesInLayer_/2;

		SizeType type = (tau_ % 3);
		SizeType offset = 0;
		if (type == 2) offset = 1;
		mu_.resize(totalUs,4); // (2,2)
		for (SizeType i = 0; i < totalUs; ++i) {
			mu_(i,0) = new TensorLeg(2*i+offset,TensorLegType::IN);
			mu_(i,1) = new TensorLeg(2*i+1+offset,TensorLegType::IN);
			mu_(i,2) = new TensorLeg(2*i+offset,TensorLegType::OUT);
			mu_(i,3) = new TensorLeg(2*i+1+offset,TensorLegType::OUT);
		}

		SizeType totalVs = totalUs;
		mw_.resize(totalVs,3); // (2,1)
		for (SizeType i = 0; i < totalVs; ++i) {
			SizeType first = 2*i + 1;
			SizeType second = 2*i + 2;
			if (type == 1 && i == 0) {
				first = sitesInLayer_;
				second = 2*i;
			}
			if (type == 1 && i > 0) {
				first = 2*i-1;
				second = 2*i;
			}
			if (type == 2 && i == 0) {
				first = sitesInLayer_;
				second = 2*i+1;
			}
			if (type == 2 && i > 0) {
				first = 2*i;
				second = 2*i+1;
			}
			mw_(i,0) = new TensorLeg(first,TensorLegType::IN);
			mw_(i,1) = new TensorLeg(second,TensorLegType::IN);
			mw_(i,2) = new TensorLeg(i,TensorLegType::OUT);
			outputSites_++;
		}
	}


	void setMeraArchitecture1dTernary()
	{
		SizeType type = (tau_ % 3);
		SizeType offset = 0;
		if (type == 0) offset = 1;
		if (type == 2) offset = 2;
		SizeType totalUs = (sitesInLayer_ + 1)/3;

		mu_.resize(totalUs,4); // (2,2) == (ins, outs)
		for (SizeType i = 0; i < totalUs; ++i) {
			SizeType first = 3*i + offset - 1;
			if (i == 0 && type == 1) first = sitesInLayer_;
			mu_(i,0) = new TensorLeg(first,TensorLegType::IN);
			mu_(i,1) = new TensorLeg(3*i + offset,TensorLegType::IN);
			mu_(i,2) = new TensorLeg(first,TensorLegType::OUT);
			mu_(i,3) = new TensorLeg(3*i + offset,TensorLegType::OUT);
		}

		SizeType totalVs = (sitesInLayer_ + 1)/3 + type - 1;
		offset = 2 - type;
		mw_.resize(totalVs,4); // (3,1) == (ins, outs)
		for (SizeType i = 0; i < totalVs; ++i) {
			if (3*i + offset >= sitesInLayer_) break;
			SizeType first = 3*i + offset - 1;
			if (i == 0 && type == 2) first = sitesInLayer_;
			mw_(i,0) = new TensorLeg(first,TensorLegType::IN);
			mw_(i,1) = new TensorLeg(3*i + offset,TensorLegType::IN);
			mw_(i,2) = new TensorLeg(3*i + offset + 1,TensorLegType::IN);
			mw_(i,3) = new TensorLeg(i,TensorLegType::OUT);
			outputSites_++;
		}
	}

	void markBogusUoutputs()
	{
		SizeType totalUs = mu_.n_row();
		SizeType cols = mu_.n_col();
		PairSizeType result;
		for (SizeType i = 0; i < totalUs; ++i) {
			for (SizeType j = 0; j < cols; ++j) {
				TensorLegType* ptr = mu_(i,j);
				if (ptr == 0) continue;
				if (ptr->inOrOut != TensorLegType::OUT) continue;
				SizeType site = ptr->site;
				if (siteIsInLegOfSomeTensor(result,mw_,site,TensorLegType::IN))
					continue;
				ptr->site = sitesInLayer_;
			}
		}
	}

	bool siteIsInLegOfSomeTensor(PairSizeType& result,
	                             const MatrixTensorLegType& m,
	                             SizeType site,
	                             TensorLegType::TensorMeraType inOrOut) const
	{
		SizeType totalVs = m.n_row();
		SizeType cols = m.n_col();
		for (SizeType i = 0; i < totalVs; ++i) {
			for (SizeType j = 0; j < cols; ++j) {
				TensorLegType* ptr = m(i,j);
				if (ptr == 0) continue;
				if (ptr->inOrOut != inOrOut) continue;
				result = PairSizeType(i,j);
				if (site == ptr->site) return true;
			}
		}

		return false;
	}

	void createSrep(PsimagLite::String& srep)
	{
		VectorSizeType sindex(sitesInLayer_,0);
		SizeType scount = 0;
		SizeType fcount = 0;
		SizeType soffset = (prevLayer_) ? prevLayer_->scount_ : 0;
		PairSizeType result;

		// set sindex
		for (SizeType i = 0; i < mw_.n_row(); ++i) {
			SizeType ins = findInsOrOuts(mw_,i,TensorLegType::IN);
			// loop over ins only here
			for (SizeType j = 0; j < ins; ++j) {
				SizeType site = mw_(i,j)->site;
				if (site >= sitesInLayer_) continue;

				if (siteIsInLegOfSomeTensor(result,mu_,site,TensorLegType::OUT)) {
					sindex[site] = scount + soffset;
					scount++;
					continue;
				}
			}
		}

		for (SizeType i = 0; i < vecForUpdateOrder_.size(); ++i) {
			SizeType s = vecForUpdateOrder_[i].first;
			TensorTypeEnum t = vecForUpdateOrder_[i].second;
			PsimagLite::String ts = (t == TENSOR_TYPE_U) ? "u" :"w";
			srep += ts + ":" + ttos(tau_) + ":" + ttos(s);
			const MatrixTensorLegType& m = (t == TENSOR_TYPE_U) ? mu_ : mw_;
			SizeType ins = findInsOrOuts(m,s,TensorLegType::IN);
			SizeType outs = findInsOrOuts(m,s,TensorLegType::OUT);
			if (t == TENSOR_TYPE_W && outs > 1)
				throw PsimagLite::RuntimeError("w has outs > 1\n");

			srep += ":" + ttos(ins) + ":" + ttos(outs);

			// loop over ins
			for (SizeType j = 0; j < ins; ++j) {
				SizeType site = m(s,j)->site;
				SizeType offset = (prevLayer_ && prevLayer_->prevLayer_) ?
				            prevLayer_->prevLayer_->scount_ : 0;

				if (site >= sitesInLayer_) {
					srep += ":d";
					continue;
				}

				if (t == TENSOR_TYPE_U) {
					if (!prevLayer_) {
						assert(tau_ == 0);
						srep += ":f" + ttos(fcount);
						fcount++;
						continue;
					}

					assert(tau_ > 0);
					if (siteIsInLegOfSomeTensor(result,prevLayer_->mw_,site,TensorLegType::OUT)) {
						srep += ":s" + ttos(result.first+offset);
						continue;
					}

					throw PsimagLite::RuntimeError("u with disconnected in is not a dummy?!\n");
				}

				if (t == TENSOR_TYPE_W) {
					if (siteIsInLegOfSomeTensor(result,mu_,site,TensorLegType::OUT)) {
						srep += ":s" + ttos(sindex[site] + outputSites_);
						continue;
					}

					if (!prevLayer_) {
						assert(tau_ == 0);
						srep += ":f" + ttos(fcount);
						fcount++;
						continue;
					}

					assert(tau_ > 0);
					if (siteIsInLegOfSomeTensor(result,prevLayer_->mw_,site,TensorLegType::OUT)) {
						srep += ":s" + ttos(result.first + offset);
						continue;
					}

					throw PsimagLite::RuntimeError("tau > 0 but W in leg not connected?!\n");
				}

				throw PsimagLite::RuntimeError("Unknown TENSOR_TYPE\n");
			}

			// loop over outs
			for (SizeType j = 0; j < outs; ++j) {
				SizeType site = m(s,ins+j)->site;

				if (site >= sitesInLayer_) {
					srep += ":d";
					continue;
				}

				if (t == TENSOR_TYPE_U) {

					if (siteIsInLegOfSomeTensor(result,mw_,site,TensorLegType::IN)) {
						srep += ":s" + ttos(sindex[site] + outputSites_);
						continue;
					}

					throw PsimagLite::RuntimeError("u with disconnected out is not a dummy?!\n");
				}

				if (t == TENSOR_TYPE_W) {
					assert(outs == 1);
					assert(j == 0);
					srep += ":s" + ttos(s+soffset);
					continue;
				}

				throw PsimagLite::RuntimeError("Unknown TENSOR_TYPE\n");
			}

			srep += ";";
			if (i > 0 && i % 5 == 0) srep += "\n";
		}

		scount_ = scount+outputSites_+soffset;
	}

	SizeType findInsOrOuts(const MatrixTensorLegType& m,
	                       SizeType ind,
	                       TensorLegType::TensorMeraType type) const
	{
		SizeType sum = 0;
		SizeType cols = m.n_col();
		for (SizeType j = 0; j < cols; ++j) {
			if (m(ind,j)->inOrOut == type) sum++;
		}

		return sum;
	}

	const ParametersForSolverType& params_;
	SizeType tau_;
	SizeType sitesInLayer_;
	MeraLayer* prevLayer_;
	SizeType outputSites_;
	SizeType scount_;
	PsimagLite::String srep_;
	std::vector<TensorType> u_; // disentagler
	//	std::vector<TensorType> w_; // isometries
	//	std::vector<TensorType> rho_; // density matrices
	std::vector<PairSizeTypeType> vecForUpdateOrder_;
	MatrixTensorLegType mu_;
	MatrixTensorLegType mw_;
}; //class

} //namespace

#endif // MERALAYER_H
