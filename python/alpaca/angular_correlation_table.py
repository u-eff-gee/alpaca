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
import matplotlib as mpl
import numpy as np


class AngularCorrelationTable:
    def __init__(self, angular_correlation):
        self.angular_correlation = angular_correlation

    def print(
        self,
        theta,
        phi,
        theta_labels=None,
        phi_labels=None,
        number_format="{:7.5f}",
        separator="&",
        endline=" \\\\",
    ):

        table = ""

        for i in range(len(theta)):
            for j in range(len(phi)):
                table += (
                    (
                        number_format.format(theta[i])
                        if theta_labels is None
                        else theta_labels[i]
                    )
                    + " {} ".format(separator)
                    + (
                        number_format.format(phi[j])
                        if phi_labels is None
                        else phi_labels[j]
                    )
                    + " {} ".format(separator)
                    + number_format.format(self.angular_correlation(theta[i], phi[j]))
                    + endline
                    + "\n"
                )

        return table
