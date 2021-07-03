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

from scipy.interpolate import interp1d


def find_indices_of_extrema(y):
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
    def __init__(self, interpolations):
        self.interpolations = interpolations

    def __call__(self, x):
        y = []
        for inter in self.interpolations:
            try:
                y.append(inter(x))
            except ValueError:
                pass
        return y
