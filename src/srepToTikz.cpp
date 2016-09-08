#include "MeraToTikz.h"
#include <iostream>
#include "Vector.h"
#include <cstdlib>

int main(int, char**)
{
	char c = '\0';
	PsimagLite::String srep;
	while (std::cin>>c) {
		srep += c;
	}

	Mera::MeraToTikz<double> obj(srep);
	std::cout<<obj;
}
