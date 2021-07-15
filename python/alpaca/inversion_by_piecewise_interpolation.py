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

from warnings import warn

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
    In the special case that

    ..math:: x_i = x_{i+1},

    this function will test the inequality

    ..math:: x_{i-1} > x_i < x_{j}

    for all

    ..math:: i+1 < j < N-1

    instead, until an element of the list is encountered that is different from :math:`x_i`.
    In other words, equal elements in the list are skipped in the search for extrema.
    Note that if

    ..math:: x_i = x_{N-1}

    for

    ..math:: i < N-1,

    then there is no extremum with an index larger or equal to :math:`i`.

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
    last = 0
    this = 1
    next = 2
    while next < len(y):
        if y[this] == y[next]:
            next += 1
            continue
        if y[this] > y[last]:
            if y[this] > y[next]:
                indices.append(this)
        elif y[this] < y[last]:
            if y[this] < y[next]:
                indices.append(this)
        last = this
        this = next
        next += 1

    return indices


def safe_interp1d(x, y, kind="cubic"):
    """Wrapper of scipy.interpolate.interp1d which automatically falls back to kind='linear'

    Given a number of points larger than 3, this function works just like interp1d if only the
    'kind' optional argument was available.
    If the number of points is less than 4, the default cubic interpolation is not possible,
    and the function falls back to a linear interpolation.
    If the number of points is less than two, the usual error message of interp1d is thrown.

    Parameters
    ----------
    x, y: (N,) array_like
        See scipy.interpolate.interp1d.
    kind: str or int
        See scipy.interpolate.interp1d.

    Warns
    -----
    UserWarning
        If the number of points to interpolate is less than 4.

    Returns
    -------
    _Interp1D object
        See scipy.interpolate.interp1d.
    """
    if len(x) < 4:
        warn(
            "Interpolation of less than 4 points requested. Falling back to linear interpolation."
        )
        return interp1d(x, y, kind="linear", bounds_error=True)
    return interp1d(x, y, kind=kind, bounds_error=True)


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
        # Note that the definition of find_indices_of_extreme ensures that:
        # * There are no extrema at the limits of the array, y[0] and y[-1]
        # * extrema[i+1] > extrema[i].
        # Therefore, any of the slices below will always yield at least a length-2 array to be used
        # by safe_interp1d, and an interpolation should always be possible.
        interpolations.append(
            safe_interp1d(
                y[0 : extrema[0] + 1],
                x[0 : extrema[0] + 1],
                kind=kind,
            )
        )
        for i in range(1, len(extrema)):
            interpolations.append(
                safe_interp1d(
                    y[extrema[i - 1] : extrema[i] + 1],
                    x[extrema[i - 1] : extrema[i] + 1],
                    kind=kind,
                )
            )
        interpolations.append(
            safe_interp1d(y[extrema[-1] :], x[extrema[-1] :], kind=kind)
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
                # When a scipy.polyint._Interpolator1D object, as returned by
                # scipy.interpolate.interp1d, is called with a scalar argument, the argument is
                # wrapped in a numpy array before evaluation.
                # This yields to the inconsistent behavior that a call with a float returns an
                # np.array(float), i.e. a zero-dimensional array.
                # Using '[()]' after evaluating the interpolation unpacks the numpy array.
                y.append(inter(x)[()])
            except ValueError:
                pass
        return y
