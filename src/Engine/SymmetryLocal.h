#ifndef SYMMETRYLOCAL_H
#define SYMMETRYLOCAL_H
#include "Matrix.h"
#include "Vector.h"
#include "IoSimple.h"
#include "Sort.h"

namespace Mera {

class SymmetryLocal {

	static const SizeType MAX_LEGS = 8;

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<VectorSizeType*>::Type VectorVectorSizeType;
	typedef PsimagLite::Matrix<VectorSizeType*> MatrixOfQnsType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	SymmetryLocal(SizeType ntensors, const VectorSizeType& qOne)
	    : qOne_(qOne), matrix_(ntensors, MAX_LEGS), nameId_(ntensors,"")
	{}

	SymmetryLocal(PsimagLite::String filename)
	{
		PsimagLite::IoSimple::In io(filename);

		io.read(qOne_,"qOne");

		int total = 0;
		io.readline(total, "SymmTensors=");
		if (total < 1)
			throw PsimagLite::RuntimeError("SymmetryLocal: reading file failed\n");

		nameId_.resize(total);
		matrix_.setTo(0);
		matrix_.resize(total, MAX_LEGS);
		for (int i = 0; i < total; ++i) {
			PsimagLite::String tmp;
			io.readline(tmp, "SymmForTensor=");
			nameId_[i] = tmp;

			int x = 0;
			io.readline(x, "Total=");
			if (x < 1)
				throw PsimagLite::RuntimeError("SymmetryLocal: reading file failed\n");
			for (int j = 0; j < x; ++j) {
				VectorSizeType* v = new VectorSizeType;
				io.read(*v, "Leg" + ttos(j));
				matrix_(i,j) = v;
				garbage_.push_back(v);
			}
		}
	}

	~SymmetryLocal()
	{
		for (SizeType i = 0; i < garbage_.size(); ++i) {
			delete garbage_[i];
			garbage_[i] = 0;
		}
	}

	void save(std::ostream& os) const
	{
		os<<"qOne\n";
		os<<vectorToString(qOne_);
		os<<"\n";
		SizeType n = matrix_.n_row();
		SizeType m = matrix_.n_col();
		os<<"SymmTensors="<<effectiveTensors()<<"\n";
		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String str("");
			SizeType count = 0;
			for (SizeType j = 0; j < m; ++j) {
				if (matrix_(i,j) == 0) continue;
				const VectorSizeType& v = *(matrix_(i,j));
				if (v.size() == 0) continue;
				str += "Leg" + ttos(j) + " ";
				str += vectorToString(v);
				str += "\n";
				++count;
			}

			if (count == 0) continue;
			os<<"SymmForTensor="<<nameId_[i]<<"\n";
			os<<"Total="<<ttos(count)<<"\n";
			os<<str;
		}
	}

	void addTensor(PsimagLite::String str,
	               VectorVectorSizeType& q,
	               const VectorSizeType& iperm)
	{
		nameId_.push_back(str);
		indexOfNameId_[str] = nameId_.size() - 1;
		SizeType nrow = matrix_.n_row();
		SizeType ncol = matrix_.n_col();
		MatrixOfQnsType m(nrow + 1, ncol);
		for (SizeType i = 0; i < nrow; ++i)
			for (SizeType j = 0; j < ncol; ++j)
				m(i, j) = matrix_(i, j);

		for (SizeType i = 0; i < q.size(); ++i)
			m(nrow, i) = q[iperm[i]];

		matrix_ = m;
	}

	void setNameId(SizeType i, PsimagLite::String str)
	{
		assert(i < nameId_.size());
		nameId_[i] = str;
	}

	void setQ(SizeType tensorIndex, SizeType legTag)
	{
		matrix_(tensorIndex,legTag) = &qOne_;
	}

	void setQ(SizeType tensorIndex, SizeType legTag, VectorSizeType* q)
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

	void addIdentity(SizeType id, SizeType dim)
	{
		VectorVectorSizeType q(2, static_cast<VectorSizeType*>(0));
		VectorSizeType* qq = new VectorSizeType(dim, 0);
		PsimagLite::String str("i");
		str += ttos(id);
		for (SizeType i = 0; i < dim; ++i)
			(*qq)[i] = i;

		q[0] = qq;
		garbage_.push_back(qq);

		VectorSizeType* qq2 = new VectorSizeType(dim, 0);
		(*qq2) = (*qq);

		q[1] = qq2;
		garbage_.push_back(qq2);

        addTensor(str, q, *qq);
	}

	SizeType size() const { return nameId_.size(); }

	VectorSizeType* q(SizeType tensorIndex, SizeType legTag)
	{
		assert(matrix_(tensorIndex,legTag));
		return matrix_(tensorIndex,legTag);
	}

	SizeType nameIdToIndex(PsimagLite::String str) const
	{
		if (indexOfNameId_.size() == 0)
			nameIdToIndex();

		if (indexOfNameId_.count(str) == 0)
			return nameId_.size();

		SizeType ind = indexOfNameId_[str];
		assert(ind < nameId_.size() && nameId_[ind] == str);
		return ind;
	}

	const VectorSizeType& qOne() const { return qOne_; }

	static SizeType truncateDimension(const VectorSizeType& dim, SizeType m)
	{
		SizeType x = productOf(dim);
		if (m == 0) return x;
		return std::min(m,x);
	}

private:

	void nameIdToIndex() const
	{
		SizeType n = nameId_.size();
		for (SizeType i = 0; i < n; ++i)
			indexOfNameId_[nameId_[i]] = i;
	}

	SizeType effectiveTensors() const
	{
		SizeType n = matrix_.n_row();
		SizeType m = matrix_.n_col();
		SizeType counter = 0;
		for (SizeType i = 0; i < n; ++i) {
			SizeType count = 0;
			for (SizeType j = 0; j < m; ++j) {
				if (matrix_(i,j) == 0) continue;
				const VectorSizeType& v = *(matrix_(i,j));
				if (v.size() == 0) continue;
				++count;
			}

			if (count == 0) continue;
			++counter;
		}

		return counter;
	}

	PsimagLite::String vectorToString(const VectorSizeType& v) const
	{
		SizeType n = v.size();
		PsimagLite::String str(ttos(n) + "   ");
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

#if 0
		for (SizeType i = 0; i < discarded; ++i) {
			r[i] = static_cast<SizeType>(drand48()*n);
			while (std::find(r.begin(), r.begin() + i, r[i]) != r.begin() + i)
				r[i] = static_cast<SizeType>(drand48()*n);
		}
#else
		VectorSizeType sorted = qq;
		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType iperm(qq.size());
		sort.sort(sorted,iperm);
		SizeType low = discarded/2;
		SizeType high = qq.size() - discarded/2;
		if (discarded & 1) ++low;
		for (SizeType i = 0; i < low; ++i)
			r[i] = iperm[i];

		for (SizeType i = high; i < qq.size(); ++i)
			r[low++] = iperm[i];
#endif

		VectorSizeType qqq(m, 0);
		SizeType counter = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (std::find(r.begin(), r.end(), i) != r.end())
				continue;
			assert(counter < qqq.size());
			qqq[counter++] = qq[i];
		}

		qq = qqq;
	}

	VectorSizeType qOne_;
	MatrixOfQnsType matrix_;
	VectorStringType nameId_;
	PsimagLite::Vector<VectorSizeType*>::Type garbage_;
	mutable std::map<PsimagLite::String, SizeType> indexOfNameId_;
}; // class SymmetryLocal
} // namespace Mera
#endif // SYMMETRYLOCAL_H
