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

import numpy as np

from .analyzing_power import AnalyzingPower
from .angular_correlation import AngularCorrelation
from .state import State
from .transition import Transition


def find_delta_brute_force(
    analyzing_power,
    analyzing_power_value,
    delta_values,
    theta,
    n_delta=int(1e3),
    atol=1e-3,
    abs_delta_max=100.0,
    return_intervals=False,
):

    arctan_delta_max = np.arctan(abs_delta_max)
    arctan_deltas = np.linspace(-arctan_delta_max, arctan_delta_max, n_delta)
    deltas = np.tan(arctan_deltas)

    delta_results = []
    delta_matches = [False] * n_delta

    for i, delta in enumerate(deltas):
        cascade_steps = []
        for j, cas_ste in enumerate(analyzing_power.angular_correlation.cascade_steps):
            if isinstance(delta_values[j], str):
                cascade_steps.append(
                    [
                        Transition(
                            cas_ste[0].em_char,
                            cas_ste[0].two_L,
                            cas_ste[0].em_charp,
                            cas_ste[0].two_Lp,
                            delta,
                        ),
                        cas_ste[1],
                    ]
                )
            else:
                cascade_steps.append(
                    [
                        Transition(
                            cas_ste[0].em_char,
                            cas_ste[0].two_L,
                            cas_ste[0].em_charp,
                            cas_ste[0].two_Lp,
                            delta_values[j],
                        ),
                        cas_ste[1],
                    ]
                )
        ana_pow = AnalyzingPower(
            AngularCorrelation(
                analyzing_power.angular_correlation.initial_state, cascade_steps
            ),
            convention=analyzing_power.convention,
        )(theta)

        if isinstance(analyzing_power_value, (int, float)):
            if np.abs(ana_pow - analyzing_power_value) < atol:
                delta_results.append(delta)
                delta_matches[i] = True
        else:
            if (
                analyzing_power_value[0] - atol
                <= ana_pow
                <= analyzing_power_value[1] + atol
            ):
                delta_results.append(delta)
                delta_matches[i] = True

    if return_intervals:
        intervals = []
        interval_start = None
        for i, matching in enumerate(delta_matches):
            if matching and interval_start is None:
                interval_start = deltas[i]
            if not matching and interval_start is not None:
                intervals.append([interval_start, deltas[i]])
                interval_start = None
            if i == n_delta - 1 and interval_start is not None:
                intervals.append([interval_start, deltas[-1]])
        return intervals

    return (delta_results, delta_matches)
