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

from alpaca.analyzing_power import AnalyzingPower
from alpaca.angular_correlation import AngularCorrelation
from alpaca.find_delta import find_delta_brute_force
from alpaca.state import POSITIVE, State
from alpaca.transition import ELECTRIC, MAGNETIC, Transition


def test_find_delta():

    delta = 0.5

    ana_pow = AnalyzingPower(
        AngularCorrelation(
            State(0, POSITIVE),
            [
                [Transition(MAGNETIC, 2, ELECTRIC, 4), State(2, POSITIVE)],
                [Transition(MAGNETIC, 2, ELECTRIC, 4, delta), State(2, POSITIVE)],
            ],
        )
    )

    theta = 0.5 * np.pi

    ana_pow_val = ana_pow(theta)

    delta_results = find_delta_brute_force(ana_pow, ana_pow_val, (0.0, "delta"), theta)

    assert len(delta_results) == 2
    assert np.isclose([delta], [delta_results[0]], atol=1e-2)
