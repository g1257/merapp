#ifndef PROGRAMGLOBALS_H
#define PROGRAMGLOBALS_H

namespace Mera {

struct ProgramGlobals {

	static SizeType packTimeSpace(SizeType time, SizeType space, SizeType tauMax)
	{
		return time + space*tauMax;
	}

	static void unpackTimeAndSpace(SizeType& time,
	                               SizeType& space,
	                               SizeType id,
	                               SizeType tauMax)
	{
		time = id % tauMax;
		space = id/tauMax;
	}

	static PsimagLite::String addLf(PsimagLite::String str, SizeType each)
	{
		SizeType l = str.length();
		PsimagLite::String str2("");
		for (SizeType i = 0; i < l; ++i) {
			if (i > 0 && i % each == 0) str2 += "\n";
			str2 += str[i];
		}

		return str2;
	}
};
} // namespace Mera
#endif // PROGRAMGLOBALS_H
