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

/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef MERA_INPUT_CHECK_H
#define MERA_INPUT_CHECK_H
#include <vector>
#include <stdexcept>
#include "Options.h"
#include "Tokenizer.h"

namespace Mera {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	InputCheck() : optsReadable_(0)
	{
		PsimagLite::String knownLabels = "TauMax IterMera IterTensor DimensionSrep";
		knownLabels += " TensorId Layer IgnoreTerm Environ Terms";
		PsimagLite::tokenizer(knownLabels,knownLabels_," ");
	}

	~InputCheck()
	{
		if (optsReadable_!=0) delete optsReadable_;
	}

	bool check(const PsimagLite::String& label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType line) const
	{
		if (label=="JMVALUES") {
			if (vec.size()!=2) return error1("JMVALUES",line);
			return true;
		} else if (label=="RAW_MATRIX" || label == "SpinOrbit") {
			if (!checkForMatrix(vec)) return error1(label,line);
			return true;
		} else if (label=="Connectors") {
			if (!checkForMatrix(vec) && !checkForVector(vec))
				return error1(label,line);
			return true;
		} else if (label == "MagneticField") {
			return true;
		} else if (label=="FiniteLoops") {
			SizeType n = atoi(vec[0].c_str());
			if (vec.size()!=3*n+1)  return error1("FiniteLoops",line);
			return true;
		}

		return false;
	}

	bool checkSimpleLabel(const PsimagLite::String& label,
	                      SizeType line) const
	{
		for (SizeType i = 0; i < knownLabels_.size(); ++i)
			if (knownLabels_[i] == label) return true;
		PsimagLite::String msg("WARNING: Unknown label " + label +"\n");
		std::cout<<msg;
		std::cerr<<msg;
		return false;
	}

	/* PSIDOC dmrgSolverOptions
	   \verb!SolverOptions=! in the input file must contain
		  a comma-separated list of strings. At least one of the following strings must
		  be provided:
		\begin{itemize}
			\item[none]  Use this when no options are given, because the list of
		   strings must be non-null.
				Note that ``none'' does not disable other options.

			 \item[useSu2Symmetry] Use the SU(2) symmetry for the model, and
			interpret quantum
				 numbers in the line ``QNS'' appropriately.

			 \item[nofiniteloops]  Don't do finite loops, even if provided under
			``FiniteLoops'' below.
			\item[restart] Restart from a previously saved run. See FIXME
			\item[debugmatrix] Print Hamiltonian matrix for targeted sector of
			superblock
		\end{itemize}
		*/
	void check(const PsimagLite::String& label,const PsimagLite::String& val,SizeType)
	{
		if (label!="SolverOptions") return;
		PsimagLite::Vector<PsimagLite::String>::Type registerOpts;

		registerOpts.push_back("restart");

		PsimagLite::Options::Writeable
		        optWriteable(registerOpts,PsimagLite::Options::Writeable::PERMISSIVE);
		optsReadable_ = new  OptionsReadableType(optWriteable,val);
	}

	bool isSet(const PsimagLite::String& thisOption) const
	{
		return optsReadable_->isSet(thisOption);
	}

	void checkForThreads(SizeType nthreads) const
	{
		if (nthreads==1) return;

		PsimagLite::String message1(__FILE__);
		message1 += " FATAL: You are requesting nthreads>0 but you ";
		message1 += "did not compile with USE_PTHREADS enabled\n";
		message1 += " Either set Threads=1 in the input file (you won't ";
		message1 += "have threads though) or\n";
		message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile ";
		message1 += "and recompile\n";
		throw PsimagLite::RuntimeError(message1.c_str());
	}

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr<<"USAGE is "<<name<<"\n";
	}

private:

	bool checkForVector(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() == 0) return false;
		SizeType n = atoi(vec[0].c_str());
		return (vec.size() == n+1);
	}

	bool checkForMatrix(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() < 2) return false;
		SizeType row = atoi(vec[0].c_str());
		SizeType col = atoi(vec[1].c_str());
		SizeType n = row*col;
		return (vec.size() == n+2);
	}

	bool error1(const PsimagLite::String& message,SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());

	}

	OptionsReadableType* optsReadable_;
	VectorStringType allowedFileOptions_;
	VectorStringType knownLabels_;
}; // class InputCheck
} // namespace Mera

/*@}*/
#endif

