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

#ifndef MERA_MODEL_SELECTOR_H
#define MERA_MODEL_SELECTOR_H
#include "../Models/Heisenberg/Heisenberg.h"
#include "../Models/Hubbard/Hubbard.h"

namespace Mera {

template<typename ModelBaseType>
class ModelSelector {

	typedef typename ModelBaseType::VectorType VectorType;

	// start models here:
	typedef Heisenberg<ModelBaseType> ModelHeisenbergType;
	typedef Hubbard<ModelBaseType> ModelHubbardType;
	// end models

public:

	ModelSelector(const PsimagLite::String& name, const VectorType& hTerms)
	    : name_(name),model_(0)
	{
		if (name_ == "Heisenberg") {
			model_ = new ModelHeisenbergType(hTerms);
		} else if (name_ == "Hubbard") {
			model_ = new ModelHubbardType(hTerms);
		} else {
			PsimagLite::String s(__FILE__);
			s += " Unknown model " + name_ + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	~ModelSelector()
	{
		if (model_) delete model_;
	}

	const ModelBaseType& operator()()
	{
		return *model_;
	}

private:

	ModelSelector(const ModelSelector&);

	ModelSelector& operator=(const ModelSelector&);

	PsimagLite::String name_;
	ModelBaseType* model_;

}; // ModelSelector

} // namespace Mera

/*@}*/
#endif // MERA_MODEL_SELECTOR_H

