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

import matplotlib.pyplot as plt

from alpaca.angular_correlation import AngularCorrelation
from alpaca.angular_correlation_plotter import AngularCorrelationPlotter
from alpaca.state import NEGATIVE, POSITIVE, State


def test_angular_correlation_plotter():
    ang_cor_plot = AngularCorrelationPlotter(
        AngularCorrelation(State(0, POSITIVE), [State(2, NEGATIVE), State(0, POSITIVE)])
    )

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ang_cor_plot.plot(ax)
    plt.savefig("test_angular_correlation_plotter.pdf")
