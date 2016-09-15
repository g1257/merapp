#include "Vector.h"
#include "TensorOptimizer.h"

int main(int argc, char** argv)
{
	typedef Mera::TensorOptimizer<double> TensorOptimizerType;
	typedef TensorOptimizerType::IoInType IoInType;
	PsimagLite::String file = "meraEnviron.txt";
	if (argc >= 2) file = argv[1];

	IoInType io(file);
	SizeType tauMax = 0;
	io.readline(tauMax,"TauMax=");

	TensorOptimizerType to(io,0);

	to.optimize();
}
