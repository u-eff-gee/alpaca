# **alpaca**: a linearly-polarized angular-correlation application

An object-oriented C++ library with python bindings for direction-direction and direction-polarization correlations of two photons in a cascade of transitions between (nuclear) quantum states.

## Description

The direction- and polarization vectors of different photons which are absorbed/emitted during a single internal electromagnetic transition sequence of a nucleus are correlated and depend on the angular-momentum ('spin') quantum numbers and parities of the involved quantum states.
Since the relations are well known from the quantum theory of light, a measurement of the direction of emission or the polarization of photons with respect to some well-defined coordinate system gives a model-independent access to excited-state properties.

This code implements analytical expressions for direction-direction and direction-polarization correlations between two photons in an arbitrarily long cascade.
More explicitly, the former answers the question:

> Assume that a photon was observed which belongs to a transition between two nuclear states. What is the probability (density) for observing another photon that is emitted in the same sequence at a given polar angle relative to the direction of the first one?

The latter takes into account polarization information as well:

> Assume that a photon was observed which belongs to a transition between two nuclear states. What is the probability (density) for observing another photon that is emitted in the same sequence at a given polar angle relative to the direction and a given azimuthal angle relative to the polarization vector of the first one?

The formalism used here was taken from a book chapter by L. C. Biedenharn [1].
In particular, the convention of Biedenharn for the multipole-mixing ratio is used, where the sign of the mixing ratio depends on the notion of initial/final states and intermediate states of a cascade.

## Building (Linux)

### Prerequisites (C++)

* [CMake](https://cmake.org/)
* C++ compiler which supports at least the C++11 standard.
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
* [doxygen](https://www.doxygen.nl/) and all its requirements for [displaying formulas](https://www.doxygen.nl/manual/formulas.html) (documentation, optional)

### Prerequisites (python)

The python bindings of this project invoke the shared C++ libraries via [ctypes](https://docs.python.org/3/library/ctypes.html#module-ctypes), which means that all mandatory prerequisites for the C++ code are also prerequisites for the python code.
In addition, the following libraries/tools are required:

* [python3](https://www.python.org/) (it is assumed that the command-line executable is called 'python3')
* [numpy](https://numpy.org/)
* [pytest](https://docs.pytest.org/) (testing, optional)
* [pytest-cov](https://pypi.org/project/pytest-cov/) (testing, optional)

### Build (C++)

In the following, it is assumed that `ALPACA_SOURCE_DIR` is the directory that contains the top-level `CMakeLists.txt` and this README.
Furthermore, it is assumed that `ALPACA_BUILD_DIR` is the directory where the user wants to build the code.

First, go to the build directory and configure the code:

```
$ cd ALPACA_BUILD_DIR
$ cmake -DBUILD_DOCUMENTATION=ON ALPACA_SOURCE_DIR
```

The code above contains the optional `-DBUILD_DOCUMENTATION=ON` argument, which is the only custom build option of `alpaca` at the moment.
When activated, it will build the `doxygen` documentation of the C++ code.

After configuring, compile the code by typing:

```
$ cmake --build .
```

### Optional Testing (C++)

To run a self test of the C++ code, type:

```
$ cd ALPACA_BUILD_DIR/test
$ ctest
```

### Build (python)

Follow the steps of the previous section.
After that, the python code can be found in `ALPACA_BUILD_DIR/python`.
For a system-wide install of the library, type:

```
$ python3 setup.py install
```

in that directory.

### Optional Testing (python)

To run a self test of the python code, type:

```
$ cd ALPACA_BUILD_DIR/python
$ pytest
```

Note that the python tests only test the python API.
A detailed test of the angular correlation formalism is performed for the C++ code only.
Using `pytest-cov`, the tests create a coverage report.

## Usage

As an example, consider an experiment in which a nucleus with a 0^+ ground state is excited by a photon beam that propagates along the positive z axis and has a linear polarization along the x axis.
The corresponding electric-dipole excitation renders the nucleus in an excited state with the quantum numbers 1^-, from which it decays to a low-lying 2^+ state via a mixed electric-dipole/magnetic-quadrupole transition with a mixing ratio of 1.
The goal is to evaluate the angular correlation at a polar angle of 90 degrees and an azimuthal angle of 0 degrees.

Note that all functions take the spin quantum number times 2 as their argument.

Both the C++ code and the python code create an `AngularCorrelation` object, which is characterized by an initial state and a sequence of transition-state pairs that represents the steps of the cascade.
States and Transitions are passed as `State` and `Transition` objects.
The `AngularCorrelation` object can be called `(theta, phi)` with a polar angle `theta` and an azimuthal angle `phi`.

### Example (C++)

```
#include <gsl/gsl_math.h>

#include <iostream>

#include "AngularCorrelation.hh"
#include "State.hh"
#include "Transition.hh"

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

### Example (python)

```
import numpy as np

from alpaca.angular_correlation import AngularCorrelation
from alpaca.state import NEGATIVE, POSITIVE, State
from alpaca.transition import ELECTRIC, MAGNETIC, Transition

ang_cor = AngularCorrelation(
    State(0, POSITIVE),
    [
        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.), State(2, NEGATIVE)],
        [Transition(ELECTRIC, 2, MAGNETIC, 4, 1.), State(4, POSITIVE)],
    ],
)

print(ang_cor(0.5*np.pi, 0.))
```

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

Copyright (C) 2021 Udo Friman-Gayer (ufg@email.unc.edu)

This program was originally developed as part of the [Geant4](https://geant4.web.cern.ch/) application [nutr](https://github.com/uga-uga/nutr), published under the GNU General Public License.

## Acknowledgements

The author would like to thank C. Iliadis for enlightening discussions about the angular correlation formalism and for help with debugging the associated modules. The author would also like to thank O. Papst for helpful discussions about the angular correlation formalism and advertise OP's angular correlation code [angcorrwat](https://github.com/op3/angcorrwat). angcorrwat is complementary to the present code in the sense that it uses the python package [sympy](https://www.sympy.org/) to obtain symbolic expressions for the angular correlations.

## References

[1] L. C. Biedenharn, 'Angular Correlations in Nuclear Spectroscopy' in F. Ajzenberg-Selove (editor), 'Nuclear Spectroscopy', Part B, Academic Press New York and London (1960)

See also `ALPACA_DIR/bibliography.bib`.