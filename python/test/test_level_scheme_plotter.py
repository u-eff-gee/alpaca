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

import matplotlib.pyplot as plt

from alpaca.level_scheme_plotter import LevelSchemePlotter
from alpaca.state import NEGATIVE, PARITY_UNKNOWN, POSITIVE, State
from alpaca.transition import ELECTRIC, EM_UNKNOWN, MAGNETIC, Transition


def test_level_scheme_plotter():
    lvl_scheme_excited = LevelSchemePlotter(
        initial_state=State(0, POSITIVE),
        cascade_steps=[
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.5), State(4, POSITIVE)],
            [Transition(EM_UNKNOWN, 4, EM_UNKNOWN, 6, 0.0), State(8, PARITY_UNKNOWN)],
        ],
        delta_labels=[r"", r"$\delta_2$", r"$\delta_3$"],
        fontsize=10,
        show_polarization=[True, False, False],
    )

    lvl_scheme_ground = LevelSchemePlotter(
        initial_state=State(0, POSITIVE),
        cascade_steps=[
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.5), State(4, POSITIVE)],
            [Transition(EM_UNKNOWN, 4, EM_UNKNOWN, 6, 0.0), State(8, PARITY_UNKNOWN)],
            [Transition(EM_UNKNOWN, 8, EM_UNKNOWN, 10, 0.0), State(0, POSITIVE)],
        ],
        delta_labels=["", r"$\delta_2$", r"$\delta_3$", ""],
        fontsize=10,
        show_polarization=[True, False, False, False],
        returns_to_initial_state=True,
    )

    fig, ax = plt.subplots(2, 1, figsize=(3, 6))
    ax[0].axis("off")
    lvl_scheme_excited.plot(ax[0])
    ax[1].axis("off")
    lvl_scheme_ground.plot(ax[1])
    plt.savefig("test_level_scheme_plotter.pdf")
