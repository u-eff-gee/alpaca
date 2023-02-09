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

from pathlib import Path

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
                [Transition(MAGNETIC, 2, ELECTRIC, 4, 0.0), State(4, POSITIVE)],
            ],
        ),
        (0.0, r"\delta_2"),
        convention="KPZ",
        returns_to_initial_state=False,
        show_polarization=[True, False],
        analyzing_power_experimental=((-0.3, 0.05), (-0.1, 0.05)),
        output_file_name="analyzing_power_0_1_2.pdf",
    ),
    AnalyzingPowerPlotter(
        AngularCorrelation(
            State(3, POSITIVE),
            [
                [Transition(MAGNETIC, 2, ELECTRIC, 4, 1.0), State(5, POSITIVE)],
                [Transition(MAGNETIC, 2, ELECTRIC, 4, 1.0), State(3, POSITIVE)],
            ],
        ),
        (r"\delta_1", r"\delta_1"),
        convention="KPZ",
        returns_to_initial_state=True,
        show_polarization=[True, False],
        analyzing_power_experimental=((-0.3, 0.05), (-0.1, 0.05)),
        output_file_name="analyzing_power_15_25_15.pdf",
    ),
]


def test_analyzing_power_single_mixing():
    for ana_pow_plot in plots:
        # if not Path(ana_pow_plot.output_file_name).exists():
        ana_pow_plot.plot()
