#ifndef MERALAYER_H
#define MERALAYER_H
#include <cassert>
#include "Vector.h"
#include "Matrix.h"
#include "TensorLeg.h"

namespace Mera {

class MeraLayer {

    typedef int TensorType;
	typedef TensorLeg TensorLegType;
	typedef PsimagLite::Matrix<TensorLegType> MatrixTensorLegType;

    enum TensorTypeEnum {TENSOR_TYPE_W,TENSOR_TYPE_U};

	enum MeraArchitectureEnum {MERA_ARCH_1D_BINARY,MERA_ARCH_1D_TERNARY,
		                       MERA_ARCH_2D_2X2,MERA_ARCH_2D_3X3};

public:

    typedef std::pair<SizeType,SizeType> PairSizeType;

    MeraLayer(SizeType tau)
        : tau_(tau),vecForUpdateOrder_(u_.size())
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

	}

	SizeType tau_;
    std::vector<TensorType> u_; // disentagler
    std::vector<TensorType> w_; // isometries
    std::vector<TensorType> rho_; // density matrices
    std::vector<PairSizeType> vecForUpdateOrder_;
	MatrixTensorLegType mu_;
	MatrixTensorLegType mv_;
}; //class

} //namespace

#endif // MERALAYER_H
