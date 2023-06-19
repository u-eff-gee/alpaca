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

import pytest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from alpaca import AngularCorrelation, Parity, State
from alpaca.angular_correlation_plotter import AngularCorrelationPlotter



def test_angular_correlation_plotter(tmp_path):
    angular_correlations = [
        ["0_1_0.pdf", State(0), [State(2), State(0)]],
    ]

    for ang_cor in angular_correlations:
        ang_cor_plot = AngularCorrelationPlotter(
            AngularCorrelation(ang_cor[1], ang_cor[2])
        )

        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        # phi_thet_psi= [0.0, 0.5 * np.pi, 0.5 * np.pi]
        ang_cor_plot.plot(ax, max_abs_value=2.3)
        # ax.set_title(
        #     r"{} $\rightarrow$ {} $\rightarrow$ {}".format(
        #         ang_cor[1].tex(parity_variable_symbol=""),
        #         ang_cor[2][0].tex(parity_variable_symbol=""),
        #         ang_cor[2][1].tex(parity_variable_symbol=""),
        #     )
        # )
        ax.tick_params(pad=-3)
        ticks = [-2, -1, 0, 1, 2]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        ax.set_xlabel(r"Beam Direction $\rightarrow$", labelpad=-2)
        if ang_cor[1].parity in (Parity.positive, Parity.negative):
            ax.set_ylabel(
                r" $\leftarrow$ Polarization Plane $\rightarrow$", labelpad=-2
            )
        plt.savefig(Path(tmp_path) / ang_cor[0])
