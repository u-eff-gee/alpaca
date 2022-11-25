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

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


class AngularCorrelationPlotter:
    def __init__(self, angular_correlation):
        self.angular_correlation = angular_correlation

    def plot(
        self, axis, PhiThetaPsi=None, n_points_per_dimension=100, max_abs_value=2.0
    ):

        theta, phi = np.meshgrid(
            np.linspace(0.0, np.pi, n_points_per_dimension),
            np.linspace(0.0, 2.0 * np.pi, n_points_per_dimension),
        )

        ang_cor = self.angular_correlation(theta, phi, PhiThetaPsi=PhiThetaPsi)

        sine_theta = np.sin(theta)
        x = ang_cor * sine_theta * np.cos(phi)
        y = ang_cor * sine_theta * np.sin(phi)
        z = ang_cor * np.cos(theta)

        color_map_max = np.max(ang_cor)
        color_map_min = np.min(ang_cor)
        color_map_norm = (ang_cor - color_map_min) / (color_map_max - color_map_min)
        color_map = mpl.cm.rainbow(color_map_norm)

        axis.set_xlim(-max_abs_value, max_abs_value)
        axis.set_ylim(-max_abs_value, max_abs_value)
        axis.set_zlim(-max_abs_value, max_abs_value)
        axis.plot_surface(
            x,
            y,
            z,
            facecolors=color_map,
            rcount=n_points_per_dimension,
            ccount=n_points_per_dimension,
        )
