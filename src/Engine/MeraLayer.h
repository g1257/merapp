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

	enum TensorTypeEnum {TENSOR_TYPE_W,TENSOR_TYPE_U};

	enum MeraArchitectureEnum {MERA_ARCH_1D_BINARY,MERA_ARCH_1D_TERNARY,
		                       MERA_ARCH_2D_2X2,MERA_ARCH_2D_3X3};

public:

	typedef TensorLegType_ TensorLegType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

	MeraLayer(const ParametersForSolverType& params, SizeType tau, SizeType sitesInLayer)
	    : params_(params),
	      tau_(tau),
	      sitesInLayer_(sitesInLayer),
	      outputSites_(0),
	      vecForUpdateOrder_(u_.size())
	{
		setUpdateOrder();
		setMeraArchitecture(MERA_ARCH_1D_TERNARY);
	}

	~MeraLayer()
	{
		for (SizeType i = 0; i < mv_.n_row(); ++i) {
			for (SizeType j = 0; j < mv_.n_col(); ++j) {
				delete mv_(i,j);
				mv_(i,j) = 0;
			}
		}
	}

	SizeType size() const
	{
		return vecForUpdateOrder_.size();
	}

	SizeType outputSites() const { return outputSites_; }

	const PairSizeType& tensorToOptimize(SizeType i) const
	{
		//std::assert(i < vecForUpdateOrder_.size());
		return vecForUpdateOrder_[i];
	}

	template<typename P>
	friend std::ostream& operator<<(std::ostream& os, const MeraLayer<P>& m)
	{
		os<<"tau="<<m.tau_<<"\n";
		os<<"sitesInThisLayer="<<m.sitesInLayer_<<"\n";
		os<<"vs="<<m.mv_.n_row()<<"\n";
		for (SizeType i = 0; i < m.mv_.n_row(); ++i) {
			os<<"isometry number "<<i<<" has sites ";
			for (SizeType j = 0; j < m.mv_.n_col(); ++j) {
				TensorLegType* ptr = m.mv_(i,j);
				if (!ptr) continue;
				PsimagLite::String val = (ptr->inOrOut == TensorLegType::OUT) ? " *" : " ";
				if (ptr->inOrOut == TensorLegType::IN && ptr->site >= m.sitesInLayer_)
					os<<val<<" o";
				else
					os<<val<<ptr->site;
			}

			os<<"\n";
		}

		os<<"-------------------\n";
		return os;
	}

private:

	void setUpdateOrder()
	{
		SizeType n = vecForUpdateOrder_.size();
		SizeType nOverTwo = n/2;

		for (SizeType i=0; i<nOverTwo; ++i) {
			vecForUpdateOrder_[2*i] = PairSizeType(i,TENSOR_TYPE_W);
			vecForUpdateOrder_[2*i+1] = PairSizeType(i,TENSOR_TYPE_U);
		}
	}

	void setMeraArchitecture(MeraArchitectureEnum what)
	{
		if (what != MERA_ARCH_1D_TERNARY)
			throw PsimagLite::RuntimeError("MERA architecture not implemented\n");

		setMeraArchitecture1dTernary();
	}

	void setMeraArchitecture1dTernary()
	{
		SizeType type = (tau_ % 3);
		SizeType totalVs = (sitesInLayer_ + 1)/3 - 1 + type;
		SizeType offset = 2 - type;
		mv_.resize(totalVs,4); // (3,1)
		for (SizeType i = 0; i < totalVs; ++i) {
			if (3*i + offset >= sitesInLayer_) break;
			SizeType first = 3*i + offset - 1;
			if (i == 0 && type == 2) first = sitesInLayer_;
			mv_(i,0) = new TensorLeg(first,TensorLeg::TensorMeraType::IN);
			mv_(i,1) = new TensorLeg(3*i + offset,TensorLeg::TensorMeraType::IN);
			mv_(i,2) = new TensorLeg(3*i + offset + 1,TensorLeg::TensorMeraType::IN);
			mv_(i,3) = new TensorLeg(i,TensorLeg::TensorMeraType::OUT);
			outputSites_++;
		}
	}

	const ParametersForSolverType& params_;
	SizeType tau_;
	SizeType sitesInLayer_;
	SizeType outputSites_;
	std::vector<TensorType> u_; // disentagler
	std::vector<TensorType> w_; // isometries
	std::vector<TensorType> rho_; // density matrices
	std::vector<PairSizeType> vecForUpdateOrder_;
	MatrixTensorLegType mu_;
	MatrixTensorLegType mv_;
}; //class

} //namespace

#endif // MERALAYER_H
