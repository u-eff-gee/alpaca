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


def invert_grid(
    x,
    fx,
    y,
    atol=1e-3,
    return_intervals=False,
):
    r"""Given a range of y values of a function y=f(x), find the corresponding range(s) for x.

    This functions inverts a one-dimensional function to find all input values that agree with
    a given range of output values.
    It expects pairs of values :math:`x_i` and :math:`y_i` for :math:`0 < i < N-1`,
    :math:`N > 0`, which are assumed to be evaluations of a function :math:`f`:

    ..math:: y_i = f(x_i)

    The function :math:`f` may be surjective, i.e. multiple values of :math:`x_i` may result in
    the same :math:`y_i`.
    By looping through all pairs of values, this function finds all :math:`x_i` that result in
    a function value in the given range :math:`\left[ y_0, y_1 \right]`.
    More precisely, this functions finds all :math:`x_i` that fulfil:

    ..math:: y_0 - \Delta y \leq \underbrace{f(x_i)}_{y_i} \leq y_1 + \Delta y

    The variable `\Delta y` is the absolute tolerance for the algorithm, which can be defined by
    the user.

    By default, this function returns the complete set of all :math:`x_i` that fulfil the
    inequality above.
    If subsequent values :math:`x_i, x_{i+1}, ..., x_{j-1}, x_j`, :math:`0 \leq i \leq j \leq N-1`
    fulfil the inequality and the grid covers the function sufficiently well, then they can be
    assumed to be part of an interval on the real axis.
    An alternative return value of this function are therefore the intervals
    :math:`\left[ x_i, x_j \right]`.

    Parameters
    ----------
    x, fx: (N,1) ndarray of float
        Values of :math:`x_i` and :math:`y_i`. Both lists must have the same length.
    y: float or [float, float]
        Value of :math:`y` or range. The two limits of the range do not have to be sorted.
    atol: float
        :math:`\Delta y`, absolute tolerance for determining the numerical equality
        (default: 0.001).
    return_intervals: bool
        Determines whether a list of :math:`x_i` (False) or a list of intervals (True)
        should be returned (default: False).

    Returns
    -------
    list of float or list of [float, float]
        Depending on the `return_intervals` setting, returns a list of all :math:`x_i`
        that match the given :math:`y` (interval) or a list of intervals of matching values.
    """

    n_x = len(x)
    if np.shape(x) != (n_x,) or np.shape(fx) != (n_x,):
        raise ValueError("Both x and fx must be (N,1) arrays [np.shape(x) == (N,)].")

    if isinstance(y, (int, float)):
        y = [y, y]
    if y[1] < y[0]:
        y = [y[1], y[0]]
    inequality = (fx >= y[0] - atol) * (fx <= y[1] + atol)
    x_results = np.extract(inequality, x)

    if return_intervals:
        intervals = []
        interval_start = None
        for i, matching in enumerate(inequality):
            if matching and interval_start is None:
                interval_start = x[i]
            if not matching and interval_start is not None:
                intervals.append([interval_start, x[i - 1]])
                interval_start = None
            if i == n_x - 1 and interval_start is not None:
                intervals.append([interval_start, x[-1]])
        return intervals

    return x_results
