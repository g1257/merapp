#include "TensorEval.h"
#include "Vector.h"

int main()
{
	PsimagLite::String str = "u0(s0)u1(s0)";

	SizeType dim0 = 5;
	typedef Mera::TensorEval<double> TensorEvalType;
	typedef TensorEvalType::TensorType TensorType;
	TensorEvalType::VectorTensorType vt(2,0);

	for (SizeType i = 0; i < vt.size(); ++i) {
		vt[i] = new TensorType(dim0,1);
		vt[i]->setToRandom();
	}

	TensorEvalType::VectorPairStringSizeType idNames;
	idNames.push_back(TensorEvalType::PairStringSizeType("u",0));
	idNames.push_back(TensorEvalType::PairStringSizeType("u",1));
	TensorEvalType tensorEval(str,vt,idNames);
	TensorEvalType::VectorSizeType freeIndices;
	std::cout<<tensorEval(freeIndices)<<"\n";

	for (SizeType i = 0; i < vt.size(); ++i) {
		delete vt[i];
		vt[i] = 0;
	}
}
