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

import matplotlib.pyplot as plt
import numpy as np

from alpaca.angular_correlation import AngularCorrelation
from alpaca.state import NEGATIVE, POSITIVE, State
from alpaca.transition import ELECTRIC, MAGNETIC, Transition

from alpaca.analyzing_power_plotter import AnalyzingPowerPlotter

plots = [
    AnalyzingPowerPlotter(
        AngularCorrelation(
            State(0, POSITIVE),
            [
                [Transition(MAGNETIC, 2, ELECTRIC, 4, 0.0), State(2, POSITIVE)],
                [Transition(MAGNETIC, 2, ELECTRIC, 4, 1.0), State(2, POSITIVE)],
            ],
        ),
        (0.0, r"\delta"),
        convention="KPZ",
        output_file_name="analyzing_power_0_1_1.pdf",
    )
]


def test_analyzing_power_single_mixing():
    for ana_pow_plot in plots:
        ana_pow_plot.plot()
