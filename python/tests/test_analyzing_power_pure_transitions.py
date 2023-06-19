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

from alpaca import AngularCorrelation, Parity, State, EMCharacter, Transition

from alpaca.analyzing_power import AnalyzingPower


def test_analyzing_power_pure_transitions(tmp_path):
    spins = [2, 4, 6]
    spin_markers = ["o", "s", "^"]
    parity_colors = ["crimson", "royalblue"]
    analyzing_powers = []

    for two_J in spins:
        em = EMCharacter.electric if two_J % 4 == 0 else EMCharacter.magnetic
        emp = (
            EMCharacter.electric if em == EMCharacter.magnetic else EMCharacter.magnetic
        )

        analyzing_powers.append(
            AnalyzingPower(
                AngularCorrelation(
                    State(0, Parity.positive),
                    [
                        [
                            Transition(em, two_J, emp, two_J + 2, 0.0),
                            State(two_J, Parity.positive),
                        ],
                        [
                            Transition(em, two_J, emp, two_J + 2, 0.0),
                            State(0, Parity.positive),
                        ],
                    ],
                ),
                convention="KPZ",
            )
        )

        em = EMCharacter.electric if two_J % 4 != 0 else EMCharacter.magnetic
        emp = (
            EMCharacter.electric if em == EMCharacter.magnetic else EMCharacter.magnetic
        )

        analyzing_powers.append(
            AnalyzingPower(
                AngularCorrelation(
                    State(0, Parity.positive),
                    [
                        [
                            Transition(em, two_J, emp, two_J + 2, 0.0),
                            State(two_J, Parity.negative),
                        ],
                        [
                            Transition(em, two_J, emp, two_J + 2, 0.0),
                            State(0, Parity.positive),
                        ],
                    ],
                ),
                convention="KPZ",
            )
        )

    _fontsize_axis_label = 13
    _xy_max = 1.5

    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))
    ax.set_xlabel(r"$A(\theta = 45^\circ)$", fontsize=_fontsize_axis_label)
    ax.set_xlim(-_xy_max, _xy_max)
    ax.set_xticks([-1, -0.5, 0.0, 0.5, 1.0])
    ax.set_ylabel(r"$A(\theta = 90^\circ)$", fontsize=_fontsize_axis_label)
    ax.set_ylim(-_xy_max, _xy_max)
    for i, ana in enumerate(analyzing_powers):
        ax.plot(
            [ana(0.25 * np.pi)],
            [ana(0.5 * np.pi)],
            spin_markers[i // 2],
            color=parity_colors[0] if i % 2 == 0 else parity_colors[1],
            markersize=8,
            # label=ana.angular_correlation.cascade_steps[0][1].tex(),
        )
    ax.legend()
    plt.tight_layout()
    plt.savefig(Path(tmp_path) / "analyzing_power_0_123_0.pdf")
