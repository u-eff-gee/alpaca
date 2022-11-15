[![Build Status](https://travis-ci.com/uga-uga/alpaca.svg?branch=master)](https://app.travis-ci.com/github/uga-uga/alpaca)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# **alpaca**: a linearly-polarized angular-correlation application

An object-oriented C++ library with python bindings for direction-direction and polarization-direction correlations of two photons in a cascade of transitions between (nuclear) quantum states.

## Table of Contents

1. [Description](#1-description)
2. [Build](#2-build)
   1. [Prerequisites (C++)](#2i-prerequisites-c)
   2. [Prerequisites (python)](#2ii-prerequisites-python)
   3. [Build (C++)](#2iii-build-c)
   4. [Optional Testing (C++)](#2iv-optional-testing-c)
   5. [Build (python)](#2v-build-python)
   6. [Optional Testing (python)](#2vi-optional-testing-python)
3. [Usage](#3-usage)
   1. [Example (C++)](#3i-example-c) 
   2. [Example (python)](#3ii-example-python) 
4. [Conventions](#4-conventions)
   1. [Coordinate System](#4i-coordinate-system)
   2. [Phases of Multipole Mixing Ratios](#4ii-phases-of-multipole-mixing-ratios)
5. [License](#5-license)
6. [Acknowledgments](#6-acknowledgments)
7. [References](#7-references)

## 1. Description

The direction- and polarization vectors of different photons, which are absorbed/emitted during a single internal electromagnetic transition sequence of a nucleus, are correlated. 
The correlation depends on the angular-momentum ('spin') quantum numbers and parities of the involved quantum states.
Since the relations are well known from the quantum theory of light, a measurement of the direction of emission or the polarization of photons with respect to some well-defined coordinate system gives model-independent access to excited-state properties.

This code implements analytical expressions for direction-direction and polarization-direction correlations between two photons in an arbitrarily long cascade.
A direction-direction correlation answers the question:

> Assume that a photon was observed which belongs to a transition between two nuclear states. What is the probability (density) for observing another photon, which is emitted in the same sequence, at a given polar angle relative to the direction of the first one?

A polarization-direction correlation takes polarization information into account as well:

> Assume that a photon was observed which belongs to a transition between two nuclear states. What is the probability (density) for observing another photon, which is emitted in the same sequence, at a given polar angle relative to the direction and a given azimuthal angle relative to the polarization vector of the first one?

This code was used to produce all numerical results in our review article on the angular-correlation formalism in the European Physical Journal A [1] (see also section 'Definition of Angles' below).
It has been tested against an equivalent code in the `R` programming language by the coauthor of Ref. [1], C. Iliadis, and against many numerical results reported in the literature.

The formalism of this code, as well as the review article Ref. [1], are based a book chapter by L. C. Biedenharn [2].
In particular, the convention of Biedenharn for the multipole-mixing ratio is used, where the sign of the mixing ratio depends on the notion of initial/final states and intermediate states of a cascade.
See also [Sec. 4.ii](#4ii-phases-of-multipole-mixing-ratios) for more information.

## 2. Build

This section describes how to build, and optionally install, the `alpaca` code.
The C++ libraries can be used without building/installing the python code.
The python bindings of this project, on the other hand, invoke the shared C++ libraries via [ctypes](https://docs.python.org/3/library/ctypes.html#module-ctypes). This means that all mandatory prerequisites for the C++ libraries are also prerequisites for the python code, and the former has to be built first.

### 2.i Prerequisites (C++)

* [CMake](https://cmake.org/) (version >= 3.16 required for installation [3])
* C++ compiler which supports at least the C++11 standard.
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
* [doxygen](https://www.doxygen.nl/) and all its requirements for [displaying formulas](https://www.doxygen.nl/manual/formulas.html) and [graphs](https://www.doxygen.nl/manual/diagrams.html) (documentation, optional)

### 2.ii Prerequisites (python)

Besides all prerequisites for the C++ libraries, the following libraries/tools are required:

* [python3](https://www.python.org/) (it is assumed that the command-line executable is called 'python3')
* [numpy](https://numpy.org/)
* [matplotlib](https://matplotlib.org/)
* [black](https://black.readthedocs.io/) (testing, optional)
* [pytest](https://docs.pytest.org/) (testing, optional)
* [pytest-cov](https://pypi.org/project/pytest-cov/) (testing, optional)
* [tox](https://tox.readthedocs.io) (testing, optional)

### 2.iii Build (C++)

In the following, it is assumed that `ALPACA_ROOT_DIR` is the directory that contains the top-level `CMakeLists.txt` and this README.
Furthermore, it is assumed that `ALPACA_BUILD_DIR` is the directory where the user wants to build the code.

First, go to the build directory and configure the code:

```
$ cd ALPACA_BUILD_DIR
$ cmake -DBUILD_DOCUMENTATION=ON -DBUILD_TESTS=ON ALPACA_ROOT_DIR
```

The code above contains the optional `-DBUILD_DOCUMENTATION=ON` and `-DBUILD_TESTS=ON` arguments, which are the custom build options of `alpaca` at the moment.
When activated, the former will build the `doxygen` documentation, and the latter will build the self tests of the C++ code.
Both default to `OFF`.

After configuring, compile the code by typing:

```
$ cmake --build .
```

To install `alpaca` in the system, type:

```
$ cmake --install .
```

The installation is not required to use the C++ libraries (the location of the libraries can be passed to the compiler manually as shown below in the C++ example) or to install the python code. See also footnote [3].

#### 2.iv Optional Testing (C++)

To run a self test of the C++ code, type:

```
$ cd ALPACA_BUILD_DIR/test
$ ctest
```

### 2.v Build (python)

Follow the steps for the C++ build in the previous section.
After that, the python code can be found in `ALPACA_BUILD_DIR/python`.
For a system-wide install of the library, type:

```
$ cd ALPACA_BUILD_DIR/python
$ python3 setup.py install
```

in that directory.

#### 2.vi Optional Testing (python)

To run a self test of the python code, type:

```
$ cd ALPACA_BUILD_DIR/python
$ tox
```

The `tox` tool will run `pytest` to perform self tests and `black` to check whether the formatting of the code corresponds to `black`'s requirements.
Using `pytest-cov`, a coverage report for the self tests will be created.

Note that the python tests only ensure that the python API works.
A detailed test of the angular correlation formalism is performed for the C++ code only.

## 3. Usage

As an example, consider an experiment in which a nucleus with a 0^+ ground state is excited by a photon beam that propagates along the positive z axis and has a linear polarization along the x axis.
The corresponding electric-dipole excitation renders the nucleus in an excited state with the quantum numbers 1^-, from which it decays to a low-lying 2^+ state via a mixed electric-dipole/magnetic-quadrupole transition with a mixing ratio of 1.
The goal is to evaluate the angular correlation at a polar angle of 90 degrees and an azimuthal angle of 0 degrees.

Note that all functions take the spin quantum number times 2 as their argument.

Both the C++ code and the python code create an `AngularCorrelation` object, which is characterized by an initial state and a sequence of transition-state pairs that represents the steps of the cascade.
States and Transitions are passed as `State` and `Transition` objects.
The `AngularCorrelation` object can be called using the arguments `(theta, phi)`, i.e. a polar angle `theta` and an azimuthal angle `phi`.

### 3.i Example (C++)

```
#include <gsl/gsl_math.h> // To be able to use pi

#include <iostream>

#include "AngularCorrelation.hh"

int main(){
    AngularCorrelation ang_cor{
        State(0, positive), // Initial state
        {
            {
                Transition(electric, 2, magnetic, 4, 0.), // Excitation 
                State(2, negative)
            },
            {
                Transition(electric, 2, magnetic, 4, 1.), // Decay 
                State(4, positive)
            },
        }
    };

    std::cout << ang_cor(0.5*M_PI, 0.) << std::endl;
}
```

In order to run this code, which is assumed to be in a file called `test.cpp`, it has to be compiled against the `alpaca` libraries.
For example, when using the GNU Compiler Collection ([GCC](https://gcc.gnu.org/)), type:

```
$ g++ test.cpp -I ALPACA_ROOT_DIR/include -L ALPACA_BUILD_DIR/source -langular_correlation -Wl,-rpath,ALPACA_BUILD_DIR/source
```

This will create an executable binary for the test.

### 3.ii Example (python)

```
import numpy as np # To be able to use pi

from alpaca.angular_correlation import AngularCorrelation
from alpaca.state import NEGATIVE, POSITIVE, State
from alpaca.transition import ELECTRIC, MAGNETIC, Transition

ang_cor = AngularCorrelation(
    State(0, POSITIVE), # Initial state
    [
        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.), State(2, NEGATIVE)], # Excitation
        [Transition(ELECTRIC, 2, MAGNETIC, 4, 1.), State(4, POSITIVE)], # Decay
    ],
)

print(ang_cor(0.5*np.pi, 0.))
```

In order to run this code, which is assumed to be in a file called `test.py`, type:

```
$ python3 test.py
```

## 4. Conventions

### 4.i Coordinate System

In contrast to Ref. [1] (see, in particular, Fig. 1 therein, and Sec. 2 'Formalism for a two-step angular correlation'), the `alpaca` code uses a right-handed coordinate system in which the beam is assumed to propagate along the `z` axis, and the azimuthal angle is defined with respect to the `x` axis.
With the definition of the angles in `alpaca`, the relations between a Cartesian coordinate system with `x`, `y` and `z`, and a spherical coordinate system with an azimuthal angle Φ and a polar angle θ, are given by:

 * x = sin(θ) cos(Φ)
 * y = sin(θ) sin(Φ)
 * z = cos(θ)

 This corresponds to the more common spherical coordinate system encountered in the literature [2,4].
 In Ref. [1], the relations are:

 * z = sin(θ) cos(Φ)
 * y = sin(θ) sin(Φ)
 * x = cos(θ)

In the aforementioned section of Ref. [1], the following statements can be found:

> A linearly polarized γ-ray beam (blue arrow) is incident along the positive **x** axis ...


> The angle Φ is between the **z** axis and the projection of the direction of the second γ ray onto the **y**-**z** plane (azimuthal angle).

In `alpaca`, they need to be modified to:

> A linearly polarized γ-ray beam (blue arrow) is incident along the positive **z** axis ...

> The angle Φ is between the **x** axis and the projection of the direction of the second γ ray onto the **x**-**y** plane (azimuthal angle).

### 4.ii Phases of Multipole Mixing Ratios

The sign of a transition matrix element is not an experimental observable.
For this reason, there is an ambiguity in the sign of the multipole mixing ratio that describes the contribution of different electromagnetic multipoles to the same transition.
Historically, different conventions have been in use, which are identified by the names of the proposing authors:

- Biedenharn (used in this work)
- Krane, Steffen, and Wheeler (KSW)
- Rose and Brink

The differences between the conventions and conversion rules are listed in our review article [1].
As already mentioned, Ref. [1] and `alpaca` use the Biedenharn convention.
While it has advantages for some nuclear resonance fluorescence (NRF) applications [1], its disadvantage is the dependence of the sign of the mixing ratio on the position of the step in the cascade.
The KSW convention, on the other hand, enforces an order in the matrix elements (initial state on the right-hand side of the matrix element, "decay matrix elements") that makes them independent of the measurement process.

Recently, a review article for the NRF method was published [5].
In contrast to our article, Zilges *et al.* use the KSW convention.
For the NRF process, where the first cascade step is an excitation and not a decay, they propose to modify the expression of KSW [6] to take into account the change of sign (see their footnote 2).
With the modified expression, the mixing ratios for both steps of an "elastic" NRF cascade have the same sign.
This seems natural, because both steps are the same transition.
Please note that - in contrast to Zilges *et al.* - the manual of the Evaluated Nuclear Structure Data File (ENSDF) [7] recommends to use the original expression of KSW with two different signs in that case.
These conflicting interpretations of the KSW conventions by the NRF community [4,5] and the ENSDF [7] are potentially confusing.

As stated in Ref. [1], it is recommended to state unambiguously in publications which convention has been used.

## 5. License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

Copyright (C) 2021, 2022 Udo Friman-Gayer (udo.friman-gayer@ess.eu)

This program was originally developed as part of the [Geant4](https://geant4.web.cern.ch/) application [nutr](https://github.com/uga-uga/nutr), published under the GNU General Public License.

## 6. Acknowledgments

The author would like to thank 

- C. Iliadis for enlightening discussions about the angular correlation formalism and for help with debugging the associated modules.
- O. Papst (OP) for helpful discussions about the angular correlation formalism and testing of the code. Furthermore, the author would like to advertise OP's angular correlation code [angcorrwat](https://github.com/op3/angcorrwat). angcorrwat is complementary to the present code in the sense that it uses the python package [sympy](https://www.sympy.org/) to obtain symbolic expressions for the angular correlations.
- N. Pietralla (NP) for helpful discussions about the mixing-ratio conventions. The writing of Ref. [5] had been in progress when our review article [1] appeared. Our use of the Biedenharn convention sparked the aforementioned discussion, and NP helped us realize the discrepancy between the nuclear resonance fluorescence literature and the ENSDF recommendations (see [Sec. 4.ii](#4ii-phases-of-multipole-mixing-ratios)). Footnote 2 in Ref. [5] addresses the issue from the point of view of Zilges *et al.*.

## 7. References

[1] C. Iliadis and U. Friman-Gayer, 'Linear polarization-direction correlations in γ-ray scattering experiments', Eur. Phys. J. A **57**, 190 (2021) (https://doi.org/10.1140/epja/s10050-021-00472-1); arXiv:2104.00228 (https://arxiv.org/abs/2104.00228)

[2] L. C. Biedenharn, 'Angular Correlations in Nuclear Spectroscopy' in F. Ajzenberg-Selove (editor), 'Nuclear Spectroscopy', Part B, Academic Press New York and London (1960)

[3] It was found that the compilation does not work with CMake versions as recent as 3.10 (default on Ubuntu 18 OS). Since the code that uses the most recent CMake features is related to the installation of the C++ libraries, it can help to comment out the last few lines in `ALPACA_DIR/CMakeLists.txt`, starting from `set(_ALPACA_HEADERS ...`. With this modification, it will not be possible to install the C++ libraries in the system's default paths. However, an installation is not required to be able to use the compiled C++ libraries or the python code.

[4] U. Kneissl, H. H. Pitz, and A. Zilges, 'Investigation of Nuclear Structure by Resonance Fluorescence Scattering', Prog. Part. Nucl. Phys. **37**, 349 (1996) (https://doi.org/10.1016/0146-6410(96)00055-5)

[5] A. Zilges, D. L. Balabanski, J. Isaak, and N. Pietralla, 'Photonuclear reactions - From basic research to applications', Prog. Part. Nucl. Phys. **122**, 103903 (2022) (https://doi.org/10.1016/j.ppnp.2021.103903)

[6] K. S. Krane, R. M. Steffen, and R. M. Wheeler, 'Directional correlations of gamma radiations emitted from nuclear states oriented by nuclear reactions or cryogenic methods', Atom. Data Nucl. Data **11**, 351-406 (1973) (https://doi.org/10.1016/S0092-640X(73)80016-6)

[7] M. J. Martin, 'Phase convention for mixing ratios in electromagnetic transitions from angular correlations and angular distributions', in 'Procedures Manual for the Evaluated Nuclear Structure Data File', ed. by M. R. Bhat., The National Nuclear Data Center, Brookhaven National Laboratory, Upton, BNL-NCS-40503 Informal Report (1987) (https://doi.org/10.2172/5262270)

See also `ALPACA_DIR/bibliography.bib`.
