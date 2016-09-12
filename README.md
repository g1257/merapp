# Preliminaries
## Disclaimer and Licensing

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
The full software license for MERA++ version 1.0.0
can be found in
file LICENSE.

## Please cite this work

MERA++ is a free and open source
multi-scale entanglement renormalization Ansatz (MERA) code
for strongly correlated electrons.
The full software license for MERA++ version 0.
can be found in
file LICENSE.
You are welcomed to use it and publish data
obtained with MERA++. If you do, please cite this
work. Explain How To Cite This Work. FIXME. TBW.

## Papers used

TBW

## Code Signature 

TBW

# Building and Running MERA++

## Required Software

* GNU C++
* PsimagLite (see below)

## Optional Software

* make or gmake (only needed to use the Makefile)
* perl (may be needed to run some auxiliary script)

## Quick Start

1. Use your distribution repository tool to install gcc with support for C++,
make, perl, and git if you don't have them.

2. Issue

    cd someDirectory/

    git clone https://github.com/g1257/PsimagLite.git

    git clone https://github.com/g1257/merapp.git

3. Compile PsimagLite

    cd PsimagLite/lib/

    make -f Makefile.sample

    cd ../../

4. Now issue

    cd merapp/src

    make

5. You can run it with TBW.

