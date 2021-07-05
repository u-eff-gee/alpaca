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

from alpaca.interval_intersections import (
    intersection_of_two_intervals,
    intersection_of_interval_with_list_of_intervals,
    intersection,
)
from alpaca.angular_correlation import AngularCorrelation
from alpaca.state import POSITIVE, State
from alpaca.transition import ELECTRIC, MAGNETIC, Transition


def test_interval_intersections():
    # Test intersection of two intervals.
    a = [0.0, 0.5]
    b = [0.6, 1.0]
    assert len(intersection_of_two_intervals(a, b)) == 0

    a = [0.6, 1.0]
    b = [0.0, 0.5]
    assert len(intersection_of_two_intervals(a, b)) == 0

    a = [0.0, 0.5]
    b = [0.4, 1.0]
    assert np.allclose(intersection_of_two_intervals(a, b), [0.4, 0.5])

    a = [0.5, 0.0]
    b = [1.0, 0.4]
    assert np.allclose(intersection_of_two_intervals(a, b), [0.4, 0.5])

    a = [0.0, 0.5]
    b = [0.3, 0.4]
    assert np.allclose(intersection_of_two_intervals(a, b), [0.3, 0.4])

    a = [0.3, 0.4]
    b = [0.2, 0.5]
    assert np.allclose(intersection_of_two_intervals(a, b), [0.3, 0.4])

    a = [0.0, 0.3]
    b = [0.3, 0.6]
    assert np.allclose(intersection_of_two_intervals(a, b), [0.3, 0.3])

    a = [0.3, 0.6]
    b = [0.2, 0.5]
    assert np.allclose(intersection_of_two_intervals(a, b), [0.3, 0.5])

    assert np.allclose(intersection_of_two_intervals(a, a), a)

    # Test intersection of an interval with a list of intervals.
    b = [[0.2, 0.4], [0.5, 0.7]]
    assert np.allclose(
        intersection_of_interval_with_list_of_intervals(a, b), [[0.3, 0.4], [0.5, 0.6]]
    )

    b = [[0.4, 0.5], [0.7, 0.8]]
    assert np.allclose(
        intersection_of_interval_with_list_of_intervals(a, b), [[0.4, 0.5]]
    )

    # Test intersection of two lists of intervals.
    a = [[0.0, 0.1], [0.3, 0.4], [0.6, 1.0]]
    b = [[0.3, 0.5], [0.8, 0.9]]
    assert np.allclose(intersection_of_interval_with_list_of_intervals(a[0], b), [])
    assert np.allclose(
        intersection_of_interval_with_list_of_intervals(a[1], b), [[0.3, 0.4]]
    )
    assert np.allclose(
        intersection_of_interval_with_list_of_intervals(a[2], b), [[0.8, 0.9]]
    )
    assert np.allclose(intersection(a, b), [[0.3, 0.4], [0.8, 0.9]])
    assert np.allclose(intersection(b, a), [[0.3, 0.4], [0.8, 0.9]])
    assert np.allclose(intersection(a, a), a)
    assert np.allclose(intersection(b, b), b)
