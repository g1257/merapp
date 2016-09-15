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
/** \ingroup MERA */
/*@{*/

/*! \file main.cpp
 *
 *  The MERA main driver
 *
 */

#include <unistd.h>
#include "MeraSolver.h"
#include <fstream>
#include "MeraToTikz.h"
#include "Version.h"

void testCase(SizeType num,
              PsimagLite::String srep,
              SizeType tauMax)
{
	std::cout<<"TEST NUMBER "<<num<<"\n";
	Mera::ParametersForSolver params(tauMax);
	Mera::MeraSolver solver(srep,params);
	solver.computeGroundState();
	std::cout<<solver;
	std::cout<<"-----------------------------------\n\n\n";

	PsimagLite::String file("meraTikzTestNumber");
	file += ttos(num);
	file += ".tex";
	std::ofstream fout(file.c_str());
	Mera::MeraToTikz<double> obj(srep,params.tauMax);
	fout<<obj;
	fout.close();
}

void oldTests(char **argv)
{
	SizeType counter = 0;
	std::cout<<argv[0]<<" version "<<MERA_VERSION<<"\n";
	PsimagLite::String str1 = "u0(f0,f1|d,s8)w0(s8,s9|s0)u3(f2,f3|s9,s10)w3(s10,s11|s1)\n";
	str1 += "u6(f4,f5|s11,s12)w6(s12,s13|s2)\n";
	str1 += "u9(f6,f7|s13,s14)w9(s14,s15|s3)u12(f8,f9|s15,s16)w12(s16,s17|s4)\n";
	str1 += "u15(f10,f11|s17,s18)w15(s18,s19|s5)u18(f12,f13|s19,s20)w18(s20,s21|s6)\n";
	str1 += "u21(f14,f15|s21,s22)w21(s22,d|s7)u1(s0,s1|s27,s28)w1(d,s27|d)\n";
	str1 += "u4(s2,s3|s29,s30)w4(s28,s29|s24)u7(s4,s5|s31,s32)w7(s30,s31|s25)\n";
	str1 += "u10(s6,s7|s33,d)w10(s32,s33|s26)u2(s24,s25|s36,s37)w2(d,s36|s34)\n";
	str1 += "u5(s26,d|s38,d)w5(s37,s38|s35)r(s34,s35)\n";
	testCase(counter++,str1,3);

	PsimagLite::String str2 = "u0(f0,f1|d,s11)w0(s11,f2,s12|s0)u3(f3,f4|s12,s13)\n";
	str2 += "w3(s13,f5,s14|s1)u6(f6,f7|s14,s15)w6(s15,f8,s16|s2)u9(f9,f10|s16,s17)\n";
	str2 += "w9(s17,f11,s18|s3)u12(f12,f13|s18,s19)w12(s19,f14,s20|s4)u15(f15,f16|s20,s21)\n";
	str2 += "w15(s21,f17,s22|s5)u18(f18,f19|s22,s23)w18(s23,f20,s24|s6)\n";
	str2 += "u21(f21,f22|s24,s25)w21(s25,f23,s26|s7)u24(f24,f25|s26,s27)w24(s27,f26,s28|s8)\n";
	str2 += "u27(f27,f28|s28,s29)w27(s29,f29,s30|s9)u30(f30,f31|s30,s31)";
	str2 += "w30(s31,f32,s32|s10)u33(f33,f34|s32,d)u1(d,s0|d,s37)w1(s37,s1,s38|d)\n";
	str2 += "u4(s2,s3|s38,s39)w4(s39,s4,s40|s34)u7(s5,s6|s40,s41)w7(s41,s7,s42|s35)\n";
	str2 += "u10(s8,s9|s42,s43)w10(s43,s10,d|s36)u2(s34,s35|s46,s47)w2(d,d,s46|s44)\n";
	str2 += "w5(s47,s36,d|s45)r(s45,s44)";
	testCase(counter++,str2,3);

	Mera::TensorSrep tensorSrep(str2);
	std::cerr<<"original " + ttos(counter) + "\n";
	tensorSrep.isValid(true);
	std::cerr<<"---------------------\n\n";
	Mera::TensorSrep tensorSrep2(tensorSrep);
	tensorSrep2.conjugate();
	testCase(counter++,tensorSrep2.sRep(),3);

	PsimagLite::String str3 = "h0(f2,f3|f0,f1)\n";
	Mera::TensorSrep tensorSrep3(str3);
	Mera::TensorSrep::VectorSizeType indicesToContract(2,0);
	indicesToContract[1] = 1;
	tensorSrep.contract(tensorSrep3,indicesToContract);
	std::cerr<<"original contracted with h" + ttos(counter) + "\n";
	tensorSrep.isValid(true);
	std::cerr<<"---------------------\n\n";
	testCase(counter++,tensorSrep.sRep(),3);

	tensorSrep.contract(tensorSrep2);
	std::cerr<<"scalar tensor" + ttos(counter) + "\n";
	tensorSrep.isValid(true);
	std::cerr<<"---------------------\n\n";
	testCase(counter++,tensorSrep.sRep(),3);

	tensorSrep2 = tensorSrep;
	tensorSrep2.eraseTensor(0);
	std::cerr<<"tensor with hole at 0, " + ttos(counter) + "\n";
	tensorSrep2.isValid(true);
	std::cerr<<"---------------------\n\n";
	testCase(counter++,tensorSrep2.sRep(),3);

	tensorSrep.eraseTensor(3);
	std::cerr<<"tensor with hole at 3, " + ttos(counter) + "\n";
	tensorSrep.isValid(true);
	testCase(counter++,tensorSrep.sRep(),3);
}

int main(int argc, char **argv)
{
	if (argc == 1) {
		std::cerr<<"USAGE: "<<argv[0]<<" tauMax\n";
		return 1;
	}

	PsimagLite::String str("");

	if (argc == 2)
		str = "u0(f0,f1|d,s0)u3(f2,f3|s1,s2)w0(s0,s1|s3)w3(s2,d|s4)r(s3,s4)\n";

	if (argc > 2) {
		char c = '\0';
		while (std::cin>>c) {
			str += c;
		}
	}

	SizeType counter = 0;
	std::cout<<argv[0]<<" version "<<MERA_VERSION<<"\n";

	testCase(counter++,str,atoi(argv[1]));
}
