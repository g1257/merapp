#ifndef SYMMETRYLOCAL_H
#define SYMMETRYLOCAL_H
#include "Matrix.h"

namespace Mera {

class SymmetryLocal {

	static const SizeType MAX_LEGS = 4;

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<const VectorSizeType*>::Type VectorVectorSizeType;
	typedef PsimagLite::Matrix<const VectorSizeType*> MatrixOfQnsType;

	SymmetryLocal(SizeType ntensors)
	    : qOne_(2,0), matrix_(ntensors, MAX_LEGS)
	{
		qOne_[1] = 1; // FIXME: PICKUP MODEL DEPENCY HERE
 	}

	~SymmetryLocal()
	{
		for (SizeType i = 0; i < garbage_.size(); ++i) {
			delete garbage_[i];
			garbage_[i] = 0;
		}
	}

	void setQ(SizeType tensorIndex, SizeType legTag)
	{
		matrix_(tensorIndex,legTag) = &qOne_;
	}

	void setQ(SizeType tensorIndex, SizeType legTag, const VectorSizeType* q)
	{
		matrix_(tensorIndex,legTag) = q;
	}

	void setQ(SizeType tensorIndex,
	          SizeType legTag,
	          const VectorVectorSizeType& q,
	          const VectorSizeType& dim,
	          SizeType m)
	{
		assert(q.size() == dim.size());
		SizeType total = productOf(dim); // untrucated
		VectorSizeType* qq = new VectorSizeType(total, 0);
		garbage_.push_back(qq);
		fillProdVector(*qq,q,dim);
		truncateVector(*qq, m);
		matrix_(tensorIndex,legTag) = qq;
	}

	const VectorSizeType* q(SizeType tensorIndex, SizeType legTag) const
	{
		assert(matrix_(tensorIndex,legTag));
		return matrix_(tensorIndex,legTag);
	}

	void print(std::ostream& os) const
	{
		SizeType n = matrix_.n_row();
		SizeType m = matrix_.n_col();
		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String str("");
			for (SizeType j = 0; j < m; ++j) {
				if (matrix_(i,j) == 0) continue;
				const VectorSizeType& v = *(matrix_(i,j));
				if (v.size() == 0) continue;
				str += "\tLeg " + ttos(j) + " indices ";
				str += vectorToString(v);
				str += "\n";
			}

			if (str == "") continue;
			os<<"Tensor "<<i<<"\n";
			os<<str;
		}
	}

	static SizeType truncateDimension(const VectorSizeType& dim, SizeType m)
	{
		SizeType x = productOf(dim);
		if (m == 0) return x;
		return std::min(m,x);
	}

private:

	PsimagLite::String vectorToString(const VectorSizeType& v) const
	{
		SizeType n = v.size();
		PsimagLite::String str("");
		for (SizeType i = 0; i < n; ++i) {
			str += ttos(v[i]) + " ";
		}

		return str;
	}

	SizeType truncateDimension(SizeType x, SizeType m) const
	{
		if (m == 0) return x;
		return std::min(m,x);
	}

	static SizeType productOf(const VectorSizeType& dim)
	{
		SizeType n = dim.size();
		if (n == 0) return 0;
		SizeType ret = dim[0];
		for (SizeType i = 1; i < n; ++i)
			ret *= dim[i];

		return ret;
	}

	void fillProdVector(VectorSizeType& qq,
	                    const VectorVectorSizeType& q,
	                    const VectorSizeType& dim) const
	{
		assert(q.size() == dim.size());
		SizeType n = qq.size();
		SizeType m = q.size();
		VectorSizeType coordinates(m,0);

		for (SizeType i = 0; i < n; ++i) {
			unpack(coordinates,i,dim);
			SizeType sum = 0;
			for (SizeType j = 0; j < m; ++j) {
				const VectorSizeType& qqq = *(q[j]);
				sum += qqq[coordinates[j]];
			}

			qq[i] = sum;
		}
	}

	void unpack(VectorSizeType& coordinates, SizeType index, const VectorSizeType& dim) const
	{
		assert(dim.size() > 0);
		SizeType tmp = index;
		SizeType n = dim.size();
		assert(coordinates.size() == n);
		for (SizeType i = 0; i < n - 1; ++i) {
			SizeType j = n - i - 1;
			div_t x = div(tmp, dim[j - 1]);
			coordinates[j] = x.quot;
			tmp = x.rem;
		}

		coordinates[0] = tmp;
	}

	void truncateVector(VectorSizeType& qq, SizeType m) const
	{
		if (m == 0) return;
		SizeType n = qq.size();
		if (n <= m) return;
		SizeType discarded = n - m;
		assert(discarded > 0);
		VectorSizeType r(discarded,0);
		for (SizeType i = 0; i < discarded; ++i) {
			r[i] = static_cast<SizeType>(drand48()*n);
			while (std::find(r.begin(), r.begin() + i, r[i]) != r.begin() + i)
				r[i] = static_cast<SizeType>(drand48()*n);
		}

		VectorSizeType qqq(m, 0);
		SizeType counter = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (std::find(r.begin(), r.end(), i) != r.end())
				continue;
			assert(counter < qqq.size());
			qqq[counter++] = i;
		}

		qq = qqq;
	}

	VectorSizeType qOne_;
	MatrixOfQnsType matrix_;
	PsimagLite::Vector<VectorSizeType*>::Type garbage_;
}; // class SymmetryLocal
} // namespace Mera
#endif // SYMMETRYLOCAL_H
