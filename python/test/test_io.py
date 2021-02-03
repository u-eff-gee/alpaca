#    This file is part of alpaca.
#
#    alpaca is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    alpaca is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with alpaca.  If not, see <https://www.gnu.org/licenses/>.
#
#    Copyright (C) 2021 Udo Friman-Gayer

import pytest

import numpy as np

from alpaca.angular_correlation import AngularCorrelation
from alpaca.state import NEGATIVE, POSITIVE, State
from alpaca.transition import ELECTRIC, MAGNETIC, Transition


def test_io():
    # The angular correlation for testing is a
    # 0^+ -> 1^- -> 0^+
    # cascade with a maximum value of 1.5 at
    # phi = pi/2, 3*pi/2
    # and a minimum value of 0. at
    # phi = 0, pi
    # for theta = pi/2.
    ang_cor = AngularCorrelation(
        State(0, POSITIVE),
        [
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(0, POSITIVE)],
        ],
    )

    theta = 0.5 * np.pi
    ang_cor_min = 0.
    ang_cor_max = 1.5
    phi_min = [0., np.pi]
    phi_max = [0.5*np.pi, 1.5*np.pi]

    # Test scalar input
    assert np.isclose(ang_cor(theta, phi_max[0]), ang_cor_max)

    # Test scalar and ndarray input
    result = ang_cor(theta, np.array([[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]]))
    assert np.allclose(result, np.array([[ang_cor_min, ang_cor_max], [ang_cor_min, ang_cor_max]]))

    result = ang_cor(np.array([[theta, theta], [theta, theta]]), phi_min[0])
    assert np.allclose(result, np.array([[ang_cor_min, ang_cor_min], [ang_cor_min, ang_cor_min]]))

    # Test ndarray input
    result = ang_cor(np.array([[theta, theta], [theta, theta]]), np.array([[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]]))
    assert np.allclose(result, np.array([[ang_cor_min, ang_cor_max], [ang_cor_min, ang_cor_max]]))

    # Error: theta has a higher dimension than phi
    with pytest.raises(ValueError):
        ang_cor(np.array([[[theta, theta], [theta, theta]]]), np.array([[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]]))
    # Error: theta has more entries than phi
    with pytest.raises(ValueError):
        ang_cor(np.array([[theta, theta], [theta, theta], [theta, theta]]), np.array([[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]]))

    # Test list input
    result = ang_cor([[theta, theta], [theta, theta]], [[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]])
    assert np.allclose(result, np.array([[ang_cor_min, ang_cor_max], [ang_cor_min, ang_cor_max]]))