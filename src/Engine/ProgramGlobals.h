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
};
} // namespace Mera
#endif // PROGRAMGLOBALS_H
