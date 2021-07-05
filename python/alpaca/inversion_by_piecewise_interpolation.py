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

import numpy as np
from scipy.interpolate import interp1d


def find_indices_of_extrema(y):
    """Find the positions of local extrema in a list of numbers

    In a list with :math:`N` entries, a local minimum is defined by

    ..math:: x_{i-1} > x_i < x_{i+1}

    ..math:: 0 < i < N-1.

    The corner points :math:`x_0` and :math:`x_{N-1}` do not count as local extrema.
    The condition for a local maximum is:

    ..math:: x_{i-1} < x_i > x_{i+1}.

    Whenever one of the two sets of inequalities is fulfilled, this function will report the
    index :math:`i`.

    Parameters:
    -----------
    y: list of int or float
        (One-dimensional) List of numbers.

    Returns:
    --------
    list of int
        List of indices of local extrema.
    """
    indices = []
    for i in range(2, len(y)):
        if y[i - 1] > y[i - 2]:
            if y[i] < y[i - 1]:
                indices.append(i - 1)
        elif y[i - 1] < y[i - 2]:
            if y[i] > y[i - 1]:
                indices.append(i - 1)

    return indices


def interpolate_and_invert(x, y, kind="cubic"):
    r"""Given pairs of values of a function, create a piecewise interpolation of its inverse

    This function inverts a one-dimensional function :math:`y = f(x)` by creating a piecewise
    interpolation :math:`\tilde{f}^{-1}` of the true inverse :math:`x = f^{-1}(y)`.
    The inverse is created from pairs of values :math:`x_i` and :math:`y_i` for
    :math:`0 < i < N-1` after finding the local extrema in the list of :math:`y` values which
    make the inverse ambiguous.
    Between every two local extrema, `scipy.interpolate.interp1d` is used to create a continuous
    mapping from :math:`y` to :math:`x` values.
    The result is return as a `PiecewiseInterpolation` object.

    Parameters
    ----------
    x, y: (N,1) ndarray of float or list of float
        Values of :math:`x_i` and :math:`y_i`. Both lists must have the same length.
    kind: string
        Type of interpolation. Equivalent to `scipy.interpolate.interp1d`'s `kind` option.
        In fact, this option is simply passed on the the `scipy` interpolator.

    Returns
    -------
    PiecewiseInterpolation
        Piecewise-interpolate approximation of the inverse of :math:`f`.

    Raises
    ------
    ValueError
        After determining the `len` of `x`, the function checks whether `x` and `fx` have the same
        shape (`numpy.shape`) `(len(x),)`
    """
    n_x = len(x)
    if np.shape(x) != (n_x,) or np.shape(y) != (n_x,):
        raise ValueError("Both x and y must be (N,1) arrays [np.shape(x) == (N,)].")

    extrema = find_indices_of_extrema(y)
    interpolations = []
    if len(extrema) == 0:
        interpolations.append(interp1d(y, x, kind=kind, bounds_error=True))
    else:
        interpolations.append(
            interp1d(
                y[0 : extrema[0]],
                x[0 : extrema[0]],
                kind=kind,
                bounds_error=True,
            )
        )
        for i in range(1, len(extrema)):
            interpolations.append(
                interp1d(
                    y[extrema[i - 1] : extrema[i]],
                    x[extrema[i - 1] : extrema[i]],
                    kind=kind,
                    bounds_error=True,
                )
            )
        interpolations.append(
            interp1d(y[extrema[-1] :], x[extrema[-1] :], kind=kind, bounds_error=True)
        )

    return PiecewiseInterpolation(interpolations)


class PiecewiseInterpolation:
    """Class for a piecewise interpolation of a function

    This class has been designed to store a list of `scipy.interpolate._Interpolator1D` objects, as
    returned by `scipy.interpolate.interp1d`.
    Interpolating a function by pieces can be helpful if the function is surjective (ambiguous, one
    x value may map to different y values).
    When the `__call__` function of a `PiecewiseInterpolation` object is used, a list of all
    possible outputs for a single input will be returned.

    It is assumed that the interpolations in the list raise a `ValueError` when an input value
    is outside their range, because `PiecewiseInterpolation` uses a `try`-`except` clause to
    check which of the pieces gives a valid result.

    Attributes
    ----------
    interpolations: list of callable objects with a single scalar input that raise a `ValueError` if the input is invalid
        Piecewise interpolations of a function.
    """

    def __init__(self, interpolations):
        """Initialize attributes

        Parameters
        ----------
        interpolations: list of callable objects with a single scalar input that raise a `ValueError` if the input is invalid
            Piecewise interpolations of a function.
        """
        self.interpolations = interpolations

    def __call__(self, x):
        """Evaluate the piecewise-interpolated function

        Parameters
        ----------
        x: float
            Input value.
        """
        y = []
        for inter in self.interpolations:
            try:
                y.append(inter(x))
            except ValueError:
                pass
        return y
