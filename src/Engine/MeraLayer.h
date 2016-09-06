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
	typedef TensorLeg TensorLegType;
	typedef PsimagLite::Matrix<TensorLegType> MatrixTensorLegType;

    enum TensorTypeEnum {TENSOR_TYPE_W,TENSOR_TYPE_U};

	enum MeraArchitectureEnum {MERA_ARCH_1D_BINARY,MERA_ARCH_1D_TERNARY,
		                       MERA_ARCH_2D_2X2,MERA_ARCH_2D_3X3};

public:

    typedef std::pair<SizeType,SizeType> PairSizeType;

    MeraLayer(const ParametersForSolverType& params)
        : params_(params),vecForUpdateOrder_(u_.size())
    {
        setUpdateOrder();
		setMeraArchitecture(MERA_ARCH_1D_TERNARY);
    }

    SizeType size() const
    {
        return vecForUpdateOrder_.size();
    }

    const PairSizeType& tensorToOptimize(SizeType i) const
    {
        //std::assert(i < vecForUpdateOrder_.size());
        return vecForUpdateOrder_[i];
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
		if (params_.tauMax & 1) setMeraArchitecture1dTernaryOdd();
		else setMeraArchitecture1dTernaryEven();
	}

	void setMeraArchitecture1dTernaryEven()
	{
		SizeType totalUs = (params_.numberOfSites + 1)/3;
		mu_.resize(totalUs,2); // only INs are accounted for for now
		for (SizeType i = 0; i < totalUs; ++i) {
			mu_(i,0) = TensorLeg(3*i,TensorLeg::TensorMeraType::IN);
			mu_(i,1) = TensorLeg(3*i + 1,TensorLeg::TensorMeraType::IN);
		}

	}

	void setMeraArchitecture1dTernaryOdd()
	{

	}

	const ParametersForSolverType& params_;
    std::vector<TensorType> u_; // disentagler
    std::vector<TensorType> w_; // isometries
    std::vector<TensorType> rho_; // density matrices
    std::vector<PairSizeType> vecForUpdateOrder_;
	MatrixTensorLegType mu_;
	MatrixTensorLegType mv_;
}; //class

} //namespace

#endif // MERALAYER_H
