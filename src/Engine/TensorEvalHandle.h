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
#ifndef TENSOREVALHANDLE_H
#define TENSOREVALHANDLE_H

namespace  Mera {

class TensorEvalHandle {

public:

	enum Status {STATUS_IDLE, STATUS_IN_PROGRESS, STATUS_DONE};

	TensorEvalHandle(Status status = STATUS_IDLE)
	    : status_(status)
	{}

	bool done() const
	{
		return (status_ == STATUS_DONE);
	}

private:

	Status status_;
};

}
#endif // TENSOREVALHANDLE_H
