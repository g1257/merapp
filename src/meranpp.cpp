#include "Vector.h"
#include "MeraSolver.h"

int main(int argc, char** argv)
{
	PsimagLite::String file = "meraEnviron.txt";
	if (argc >= 2) file = argv[1];

	Mera::MeraSolver<double> meraSolver(file);

//	TensorOptimizerType to(io,tokens[0],std::atoi(tokens[1].c_str()));

//	to.optimize(10);
}
