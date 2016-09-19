#include "Vector.h"
#include "MeraSolver.h"
#include "Version.h"
int main(int argc, char** argv)
{
	PsimagLite::String file = "meraEnviron.txt";
	if (argc >= 2) file = argv[1];

	std::cout<<"#MERA_VERSION="<<MERA_VERSION<<"\n";
	Mera::MeraSolver<double> meraSolver(file);

	meraSolver.optimize(10);
}
