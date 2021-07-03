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

# Copyright (C) 2021 Udo Friman-Gayer

import pytest

import numpy as np

from alpaca.inversion_by_piecewise_interpolation import (
    find_indices_of_extrema,
    interpolate_and_invert,
)


def f(x):
    x_squared = x * x
    return x_squared * (x_squared - 1.0)


def test_extremum_finder():

    # Test the extremum-finding method.

    extrema = find_indices_of_extrema(np.array([0.0, 1.0, 2.0, 1.0, 0.0, 1.0, -2.0]))
    assert len(extrema) == 3
    assert extrema[0] == 2
    assert extrema[1] == 4
    assert extrema[2] == 5


def test_inversion():
    x = np.linspace(-3.0, 3.0, 501)
    y = f(x)
    extrema = find_indices_of_extrema(y)
    assert len(extrema) == 3

    inverse_f = interpolate_and_invert(x, y)

    assert len(inverse_f.interpolations) == 4
    assert np.allclose(inverse_f(12.0), [-2.0, 2.0])
    assert np.allclose(inverse_f(0.0), [-1.0, 0.0, 1.0])
    assert len(inverse_f(-1.0)) == 0

    y = x
    inverse_f = interpolate_and_invert(x, y)

    assert np.isclose(inverse_f(y[0]), x[0])
    assert np.isclose(inverse_f(y[-1]), x[-1])
