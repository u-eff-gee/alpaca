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

import numpy as np


def intersection_of_two_intervals(interval_1, interval_2):
    r"""Find the intersection of two intervals

    Given two intervals :math:`\left[ a, b \right]` and :math:`\left[ c, d \right]`, this
    function finds their intersection.

    When searching for the intersection, all comparisons are made with the :math:`\leq` and
    :math:`\geq` operators.
    This means that if both arrays only have a single real number in common, i.e.
    :math:`\left[ a, b \right]` and :math:`\left[ b, d \right]`, the returned intersection will
    be :math:`\left[ b, b \right]`, and not an empty list.

    The intervals are sorted before they are processed by this function, so the interval limits can
    be unsorted.

    Parameters
    ----------
    interval_1: list or ndarray of two float
        First interval. May be unsorted.
    interval_2: list or ndarray of two float
        Second interval. May be unsorted.

    Returns
    -------
    list of two float or empty list
        Intersection of the two intervals as a new interval, or an empty list.

    Examples
    --------
    >>> intersection_of_two_intervals([0.0, 1.0], [0.5, 1.5])
    [0.5, 1.0]
    >>> intersection_of_two_intervals([0.0, 1.0], [0.5, 0.6])
    [0.5, 0.6]
    >>> intersection_of_two_intervals([0.0, 1.0], [1.0, 1.5])
    [1.0, 1.0]
    >>> intersection_of_two_intervals([0.0, 1.0], [1.5, 2.0])
    []
    """
    intersection = []

    interval_1 = np.sort(interval_1)
    interval_2 = np.sort(interval_2)

    if interval_1[0] > interval_2[1]:
        return []
    elif interval_1[0] >= interval_2[0] and interval_1[0] <= interval_2[1]:
        if interval_1[1] <= interval_2[1]:
            return [interval_1[0], interval_1[1]]
        else:
            return [interval_1[0], interval_2[1]]
    elif interval_1[1] >= interval_2[0]:
        if interval_1[1] <= interval_2[1]:
            return [interval_2[0], interval_1[1]]
        else:
            return [interval_2[0], interval_2[1]]

    return []


def intersection_of_interval_with_list_of_intervals(interval_1, list_of_intervals):
    r"""Find the intersection of an interval with a list of intervals

    Given an interval :math:`\left[ a, b \right]` and a set of intervals :math:`\left\{ \left[ c, d \right], \left[ e, f \right], ... \right \}`, this functions determines the intersections of the former with all elements of the latter.

    See also `alpaca.analyzing_power.intersection_of_two_intervals`.

    Parameters
    ----------
    interval_1: list or ndarray of two float
        First interval. May be unsorted.
    list_of_intervals: list of lists or ndarrays of two float
        List of intervals. The single intervals may be unsorted.

    Returns
    -------
    list of lists of two float, or empty list
        Intersection of the interval with the list of intervals as a new list of intervals, or an empty list.

    Examples
    --------
    >>> intersection_of_interval_with_list_of_intervals([0.0, 1.0], [[0.0, 0.5], [0.8, 1.5]])
    [[0.0, 0.5], [0.8, 1.0]]
    """
    intersections = []

    for interval in list_of_intervals:
        intersection = intersection_of_two_intervals(interval_1, interval)
        if len(intersection) > 0:
            intersections.append(intersection)

    return intersections


def intersection(list_of_intervals_1, list_of_intervals_2):
    r"""Find the intersections of a list of intervals with another list of intervals

    Given a set of intervals :math:`\left\{ \left[ \alpha, \beta \right], \left[ \gamma, \delta \right], ... \right \}` and a set of intervals :math:`\left\{ \left[ a, b \right], \left[ c, d \right], ... \right \}`, this functions determines the intersections of all of the former intervals with all elements of the latter.

    See also `alpaca.analyzing_power.intersection_of_two_intervals` and `alpaca.analyzing_power.intersection_of_interval_with_list_of_intervals`.

    Parameters
    ----------
    list_of_intervals_1: list of lists or ndarrays of two float
        First list of intervals. The single intervals may be unsorted.
    list_of_intervals_2: list of lists or ndarrays of two float
        Second list of intervals. The single intervals may be unsorted.

    Returns
    -------
    list of lists of two float, or empty list
        Intersection of the first list of intervals with the second list of intervals as a new list of intervals, or an empty list.

    Examples
    --------
    >>> intersection([[0.0, 0.1], [0.3, 0.4], [0.6, 1.0]], [[0.3, 0.5], [0.8, 0.9]])
    [[0.3, 0.4], [0.8, 0.9]]
    """
    intersections = []

    for interval_1 in list_of_intervals_1:
        intersection = intersection_of_interval_with_list_of_intervals(
            interval_1, list_of_intervals_2
        )
        if len(intersection) > 0:
            for interval in intersection:
                intersections.append(interval)

    return intersections
