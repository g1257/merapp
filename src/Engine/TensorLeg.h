#ifndef TENSORLEG_H
#define TENSORLEG_H
#include "Vector.h"

namespace Mera {

struct TensorLeg {

	enum TensorMeraType {IN, OUT};

	TensorLeg(SizeType site_, TensorMeraType inOrOut_)
	    : site(site_),inOrOut(inOrOut_)
	{}

	SizeType site;
	TensorMeraType inOrOut;
};

} //namespace Mera
#endif // TENSORLEG_H
