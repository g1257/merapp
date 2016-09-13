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
#ifndef MERALAYER_H
#define MERALAYER_H
#include <cassert>
#include "Vector.h"
#include "Matrix.h"
#include "TensorSrep.h"
#include "ProgramGlobals.h"

namespace Mera {

template<typename ParametersForSolverType>
class MeraLayer {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	enum TensorTypeEnum {TENSOR_TYPE_W,TENSOR_TYPE_U};

public:

	typedef std::pair<SizeType,TensorTypeEnum> PairSizeTypeType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

	MeraLayer(const ParametersForSolverType& params,
	          SizeType tau,
	          const TensorSrep& tensorSrep,
	          MeraLayer* prevLayer)
	    : params_(params),
	      tau_(tau),
	      tensorSrep_(tensorSrep),
	      prevLayer_(prevLayer)
	{
		setUpdateOrder(tau);
	}

	SizeType size() const
	{
		return vecForUpdateOrder_.size();
	}

	const PairSizeTypeType& tensorToOptimize(SizeType i) const
	{
		assert(i < vecForUpdateOrder_.size());
		return vecForUpdateOrder_[i];
	}

	template<typename P>
	friend std::ostream& operator<<(std::ostream& os, const MeraLayer<P>& m)
	{
		os<<"tau="<<m.tau_<<"\n";
		os<<"-------------------\n";
		return os;
	}

private:

	MeraLayer(const MeraLayer&);

	MeraLayer& operator=(const MeraLayer&);

	void setUpdateOrder(SizeType tau)
	{
		VectorSizeType mu;
		VectorSizeType mw;
		setTensorsForThisLayer(mu,mw,tau);
		SizeType n = std::min(mu.size(), mw.size());
		vecForUpdateOrder_.resize(mu.size() + mw.size());

		for (SizeType i=0; i<n; ++i) {
			vecForUpdateOrder_[2*i] = PairSizeTypeType(mu[i],TENSOR_TYPE_U);
			vecForUpdateOrder_[2*i+1] = PairSizeTypeType(mw[i],TENSOR_TYPE_W);
		}

		int x = mu.size() - mw.size();
		if (x == 0) return;

		SizeType m = std::abs(x);
		TensorTypeEnum type = (x < 0) ? TENSOR_TYPE_W : TENSOR_TYPE_U;
		const VectorSizeType& mm = (x < 0) ? mw : mu;
		for (SizeType i=n; i<n+m; ++i)
			vecForUpdateOrder_[n+i] =  PairSizeTypeType(mm[i],type);
	}

	void setTensorsForThisLayer(VectorSizeType& mu,
	                            VectorSizeType& mw,
	                            SizeType tau) const
	{
		for (SizeType i = 0; i < tensorSrep_.size(); ++i) {
			SizeType tensorX = 0;
			SizeType tensorY = 0;
			ProgramGlobals::unpackTimeAndSpace(tensorY,
			                                   tensorX,
			                                   tensorSrep_(i).id(),
			                                   params_.tauMax);
			if (tensorY != tau) continue;
			TensorStanza::TensorTypeEnum type = tensorSrep_(i).type();
			if (type == TensorStanza::TENSOR_TYPE_U)
				mu.push_back(i);
			if (type == TensorStanza::TENSOR_TYPE_W)
				mw.push_back(i);
		}
	}

	const ParametersForSolverType& params_;
	SizeType tau_;
	const TensorSrep& tensorSrep_;
	MeraLayer* prevLayer_;
	// std::vector<TensorType> u_; // disentagler
	//	std::vector<TensorType> w_; // isometries
	//	std::vector<TensorType> rho_; // density matrices
	std::vector<PairSizeTypeType> vecForUpdateOrder_;
}; //class

} //namespace

#endif // MERALAYER_H
