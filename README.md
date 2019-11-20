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

## Mission Statement

MERA++ is a C++ native application implementing the MERA algorithm for strongly correlated electron models. MERA++ is composed of two parts:
a generic engine that can handle multiple dimensions, arities, models, and geometries; and built-in models, MERA builders for different
dimensions and arities, and geometries. Features and options are chosen from a user-friendly input file.
The implementation aims to be as fast as possible, so that it compiles natively, uses optimized tensor contractions, symmetries,
layering or caching acceleration techniques, and parallelization.

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

```BASH
    $ cd someDirectory/
    $ git clone https://github.com/g1257/PsimagLite.git
    $ git clone https://github.com/g1257/merapp.git
```

3. Compile PsimagLite

```BASH
    $ cd PsimagLite
    $ git checkout features
    $ git pull origin features
    $ cd lib/
    $ ./configure.pl
    $ make -j something
    $ cd ../../
```

4. Now issue

```BASH
    $ cd merapp
    $ git checkout features
    $ git pull origin features
    $ cd src
```

Create a file myconfig.psiTag with the following.
If NOT using exatn write
```BASH
dependency exatn = CPPFLAGS += -DNO_EXATN
```
else write
```BASH
location exatn = (EXATN=/path/to/your/.exatn)
```
and then
```BASH
    $ ./configure.pl -c myconfig.psiTag
    $ make clean
    $ make -j something
```

5. You can run it with

```BASH
    $ ./merapp -n 8 -m 6  -a 2 -d 1  > out8sites
    $ ./meranpp -f out8sites.txt
```

    or run just the evaluator example driver
```BASH
    $ ./tensorEval
```

    or run from the TestSuite
```BASH
    ./meranpp -f ../TestSuite/inputs/meraEnviron1.txt
```
    etc.


