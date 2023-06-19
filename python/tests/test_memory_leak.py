# This file is part of alpaca.

# alpaca is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# alpaca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with alpaca.  If not, see <https://www.gnu.org/licenses/>.

# Copyright (C) 2021-2023 Udo Friman-Gayer

import pytest

import numpy as np

from alpaca import AngularCorrelation, State, Parity, Transition, EMCharacter
from alpaca.inversion_by_piecewise_interpolation import interpolate_and_invert
from alpaca.analyzing_power import AnalyzingPower, arctan_grid

# This test, which implements a common application of finding multipole mixing ratios that
# result in a given experimental asymmetry, was created to test for memory leaks.
# In previous versions of alpaca, the AngularCorrelation C++ object that is owned by the
# AnalyzingPower python class had not been cleaned up.
# This caused especially AnalyzingPower.evaluate() (used below) to quickly fill the available
# virtual storage, because it creates many instances of AngularCorrelation.
def test_memory_leak():
    n_monte_carlo = int(10)
    delta = arctan_grid(101)

    delta_results = []

    for _ in range(n_monte_carlo):
        ana_pow = AnalyzingPower(
            AngularCorrelation(
                State(0, Parity.positive),
                [
                    [
                        Transition(
                            EMCharacter.electric, 2, EMCharacter.magnetic, 4, 0.0
                        ),
                        State(2, Parity.negative),
                    ],
                    [
                        Transition(
                            EMCharacter.electric, 2, EMCharacter.magnetic, 4, 0.0
                        ),
                        State(4, Parity.positive),
                    ],
                ],
            )
        ).evaluate(delta, [0.0, "delta"])

        delta_inv = interpolate_and_invert(delta, ana_pow)(np.random.normal(0.0, 1.0))

        for d in delta_inv:
            delta_results.append(d)
