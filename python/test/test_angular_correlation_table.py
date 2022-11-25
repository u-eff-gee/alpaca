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

# Copyright (C) 2021, 2022 Udo Friman-Gayer

import pytest

import numpy as np

from alpaca.angular_correlation import AngularCorrelation
from alpaca.angular_correlation_table import AngularCorrelationTable
from alpaca.state import NEGATIVE, POSITIVE, State

angular_correlations = [
    ["0_1_0.txt", State(0), [State(2), State(0)]],
    ["0p_1p_0p.txt", State(0, POSITIVE), [State(2, POSITIVE), State(0, POSITIVE)]],
    ["0p_1m_0p.txt", State(0, POSITIVE), [State(2, NEGATIVE), State(0, POSITIVE)]],
    ["0_2_0.txt", State(0), [State(4), State(0)]],
    ["0_2_2.txt", State(0), [State(4), State(4)]],
]


def test_angular_correlation_table():

    theta_labels = range(0, 190, 30)
    theta = np.array(theta_labels) / 180.0 * np.pi
    theta_labels = [str(i) for i in theta_labels]
    phi_labels = range(0, 380, 30)
    phi = np.array(phi_labels) / 180.0 * np.pi
    phi_labels = [str(i) for i in phi_labels]

    for ang_cor in angular_correlations:
        ang_cor_tab = AngularCorrelationTable(
            AngularCorrelation(ang_cor[1], ang_cor[2])
        )

        with open(ang_cor[0], "w") as output_file:
            output_file.write(ang_cor_tab.print(theta, phi, theta_labels, phi_labels))
            output_file.close()
