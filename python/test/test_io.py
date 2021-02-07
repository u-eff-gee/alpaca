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
from alpaca.state import NEGATIVE, PARITY_UNKNOWN, POSITIVE, State
from alpaca.transition import ELECTRIC, EM_UNKNOWN, MAGNETIC, Transition


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
    ang_cor_min = 0.0
    ang_cor_max = 1.5
    phi_min = [0.0, np.pi]
    phi_max = [0.5 * np.pi, 1.5 * np.pi]

    # Test scalar input
    assert np.isclose(ang_cor(theta, phi_max[0]), ang_cor_max)

    # Test scalar and ndarray input
    result = ang_cor(
        theta, np.array([[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]])
    )
    assert np.allclose(
        result, np.array([[ang_cor_min, ang_cor_max], [ang_cor_min, ang_cor_max]])
    )

    result = ang_cor(np.array([[theta, theta], [theta, theta]]), phi_min[0])
    assert np.allclose(
        result, np.array([[ang_cor_min, ang_cor_min], [ang_cor_min, ang_cor_min]])
    )

    # Test ndarray input
    result = ang_cor(
        np.array([[theta, theta], [theta, theta]]),
        np.array([[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]]),
    )
    assert np.allclose(
        result, np.array([[ang_cor_min, ang_cor_max], [ang_cor_min, ang_cor_max]])
    )

    # Error: theta has a higher dimension than phi
    with pytest.raises(ValueError):
        ang_cor(
            np.array([[[theta, theta], [theta, theta]]]),
            np.array([[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]]),
        )
    # Error: theta has more entries than phi
    with pytest.raises(ValueError):
        ang_cor(
            np.array([[theta, theta], [theta, theta], [theta, theta]]),
            np.array([[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]]),
        )

    # Test list input
    result = ang_cor(
        [[theta, theta], [theta, theta]],
        [[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]],
    )
    assert np.allclose(
        result, np.array([[ang_cor_min, ang_cor_max], [ang_cor_min, ang_cor_max]])
    )

    # Test transition inference
    ang_cor = AngularCorrelation(
        State(0, POSITIVE),
        [State(2, NEGATIVE), State(0, NEGATIVE)],
    )

    assert ang_cor.cascade_steps[0][0].two_L == 2
    assert ang_cor.cascade_steps[0][0].em_char == ELECTRIC
    assert ang_cor.cascade_steps[0][0].two_Lp == 4
    assert ang_cor.cascade_steps[0][0].em_charp == MAGNETIC
    assert ang_cor.cascade_steps[0][0].delta == 0.0
    assert ang_cor.cascade_steps[1][0].two_L == 2
    assert ang_cor.cascade_steps[1][0].em_char == MAGNETIC
    assert ang_cor.cascade_steps[1][0].two_Lp == 4
    assert ang_cor.cascade_steps[1][0].em_charp == ELECTRIC
    assert ang_cor.cascade_steps[0][0].delta == 0.0

    result = ang_cor(
        [[theta, theta], [theta, theta]],
        [[phi_min[0], phi_max[0]], [phi_min[1], phi_max[1]]],
    )
    assert np.allclose(
        result, np.array([[ang_cor_min, ang_cor_max], [ang_cor_min, ang_cor_max]])
    )

    ang_cor = AngularCorrelation(
        State(0, PARITY_UNKNOWN),
        [State(2, NEGATIVE), State(0, NEGATIVE)],
    )

    assert ang_cor.cascade_steps[0][0].two_L == 2
    assert ang_cor.cascade_steps[0][0].em_char == EM_UNKNOWN
    assert ang_cor.cascade_steps[0][0].two_Lp == 4
    assert ang_cor.cascade_steps[0][0].em_charp == EM_UNKNOWN
    assert ang_cor.cascade_steps[0][0].delta == 0.0
    assert ang_cor.cascade_steps[1][0].two_L == 2
    assert ang_cor.cascade_steps[1][0].em_char == MAGNETIC
    assert ang_cor.cascade_steps[1][0].two_Lp == 4
    assert ang_cor.cascade_steps[1][0].em_charp == ELECTRIC
    assert ang_cor.cascade_steps[0][0].delta == 0.0
