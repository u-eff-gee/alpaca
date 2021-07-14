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
    safe_interp1d,
)

# This quartic function has a multiplicity-2 root at x=0, and two multiplicity-1 roots at x=-1 and
# x=1.
# There is a local maximum at x=0 with f(0) = 0 and two local minima at x=-1/sqrt(2) and
# x=1/sqrt(1) with f(+-1/sqrt(2))=-0.25.
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
    # Between -3 and 3 (actually between -infinity and infinity), the function will be split up
    # into four pieces at the three extrema.
    x = np.linspace(-3.0, 3.0, 501)
    y = f(x)
    extrema = find_indices_of_extrema(y)
    assert len(extrema) == 3

    # x and y must have the same length.
    with pytest.raises(ValueError):
        interpolate_and_invert(x, [0, 1])

    inverse_f = interpolate_and_invert(x, y)

    assert len(inverse_f.interpolations) == 4
    assert np.allclose(inverse_f(12.0), [-2.0, 2.0])
    assert np.allclose(inverse_f(0.0), [-1.0, 0.0, 1.0])
    assert len(inverse_f(-1.0)) == 0
    # The following test ensures that a scalar value is returned when the interpolated function is
    # called with a scalar argument.
    # By default, the interpolation object returned by scipy.interpolate.interp1d returns a zero-
    # dimensional numpy array when called with a scalar, which was found inconsistent, not only
    # by the author of the present code, but there are also various discussions on the web about
    # this.
    # See also alpaca/inversion_by_piecewise_interpolation.py.
    assert np.isclose(inverse_f(12.0)[0], -2.0)

    # Using the identity function, test whether the limits of the interval are treated correctly.
    y = x
    inverse_f = interpolate_and_invert(x, y)

    assert np.isclose(inverse_f(y[0]), x[0])
    assert np.isclose(inverse_f(y[-1]), x[-1])


@pytest.mark.parametrize("fx", [(range(4)), (range(3)), (range(2))])
def test_interpolation(fx):
    # Test special cases where only few points are available for the interpolation.
    # In this case, the interpolator uses a linear interpolation function, no matter what value
    # was used for the 'kind' option.
    if len(fx) < 4:
        with pytest.warns(UserWarning):
            f = safe_interp1d(range(len(fx)), fx)
    else:
        f = safe_interp1d(range(len(fx)), fx)
    assert np.isclose(f(0.5), 0.5)
