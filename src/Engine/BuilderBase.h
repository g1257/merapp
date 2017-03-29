/*
Copyright (c) 2016-2017, UT-Battelle, LLC

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
#ifndef MERA_BUILDERBASE_H
#define MERA_BUILDERBASE_H
#include "AllocatorCpu.h"
#include "TensorSrep.h"

namespace Mera {
class BuilderBase {

public:

	virtual ~BuilderBase() { }

	virtual const PsimagLite::String& srep() const = 0;

	virtual TensorSrep* buildEnergyTerm(SizeType site,
	                                    SizeType sites,
	                                    const TensorSrep& tensorSrep) const = 0;
};
}

#endif // MERA_BUILDERBASE_H
