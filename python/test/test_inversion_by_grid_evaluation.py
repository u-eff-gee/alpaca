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

from alpaca.inversion_by_grid_evaluation import invert_grid


def linear(x):
    return x


def quadratic(x):
    return x * x


def test_inversion_by_grid_evaluation():
    x = np.linspace(0.0, 1.0, 101)
    fx = linear(x)

    with pytest.raises(ValueError):
        invert_grid(x, [0, 1], 0.5)

    x_results = invert_grid(x, fx, 0.5)
    assert len(x_results) == 1
    assert x_results[0] == 0.5

    x_results = invert_grid(x, fx, 0.5, atol=0.01)
    assert len(x_results) == 3

    x_results = invert_grid(x, fx, [0.3, 0.5])
    assert len(x_results) == 21
    x_intervals = invert_grid(x, fx, [0.3, 0.5], return_intervals=True)
    assert len(x_intervals) == 1
    assert np.allclose(x_intervals, [0.3, 0.5])

    x_results = invert_grid(x, fx, [0.9, 1.1])
    assert len(x_results) == 11
    x_intervals = invert_grid(x, fx, [0.9, 1.1], return_intervals=True)
    assert len(x_intervals) == 1
    assert np.allclose(x_intervals, [0.9, 1.0])

    x_intervals = invert_grid(x, fx, [-0.1, 0.1], return_intervals=True)
    assert len(x_intervals) == 1
    assert np.allclose(x_intervals, [0.0, 0.1])

    x = np.linspace(-3.0, 3.0, 601)
    fx = quadratic(x)
    x_intervals = invert_grid(x, fx, [-1.0, 1.0], return_intervals=True)
    assert len(x_intervals) == 1
    assert np.allclose(x_intervals, [-1.0, 1.0])

    x_intervals = invert_grid(x, fx, [4.0, 1.0], return_intervals=True)
    assert len(x_intervals) == 2
    assert np.allclose(x_intervals, [[-2.0, -1.0], [1.0, 2.0]])
