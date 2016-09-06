#ifndef TENSORLEG_H
#define TENSORLEG_H
#include "Vector.h"

namespace Mera {

struct TensorLeg {

	enum TensorMeraType {IN, OUT};

	SizeType site;
	TensorMeraType inOrOut;
};

} //namespace Mera
#endif // TENSORLEG_H
