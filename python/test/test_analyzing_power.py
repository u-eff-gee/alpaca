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
from alpaca.state import NEGATIVE, POSITIVE, State
from alpaca.transition import ELECTRIC, MAGNETIC, Transition
from alpaca.analyzing_power import AnalyzingPower, arctan_grid


def test_analyzing_power():

    ang_cor = AngularCorrelation(
        State(0, POSITIVE),
        [
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(0, POSITIVE)],
        ],
    )

    ana_pow = AnalyzingPower(ang_cor)

    assert np.isclose(ana_pow(0.5 * np.pi), -1.0)

    with pytest.raises(ValueError):
        AnalyzingPower(ang_cor, PQ=2.0)

    ana_pow = AnalyzingPower(ang_cor, PQ=-0.9)

    assert np.isclose(ana_pow(0.5 * np.pi), 0.9)

    ana_pow = AnalyzingPower(ang_cor, convention="KPZ")

    assert np.isclose(ana_pow(0.5 * np.pi), 1.0)

    # Test the arbitrary-angle input

    assert np.isclose(
        ana_pow(0.1, 0.2, 0.3, 0.4),
        (ang_cor(0.2, 0.4) - ang_cor(0.1, 0.3))
        / (ang_cor(0.1, 0.3) + ang_cor(0.2, 0.4)),
    )

    # Test the arctan_grid function which creates an equidistant grid of arctan(delta) between two
    # limits.
    # In particular, test the warnings issued by this function.

    with pytest.warns(UserWarning):
        grid = arctan_grid(1)
    assert len(grid) == 2
    with pytest.warns(UserWarning):
        grid = arctan_grid(4)
    assert np.allclose(
        grid, np.tan(np.linspace(np.arctan(-100.0), np.arctan(100.0), 4))
    )

    # Test AnalyzingPower.evaluate for scalar input
    theta = 0.5 * np.pi
    ang_cor = AngularCorrelation(
        State(0, POSITIVE),
        [
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.5), State(2, NEGATIVE)],
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.5), State(0, POSITIVE)],
        ],
    )
    ana_pow = AnalyzingPower(ang_cor)

    assert ana_pow.evaluate(0.5, ["delta", "delta"], theta) == ana_pow(theta)

    # Test AnalyzingPower.evaluate for scalar input when the relation between mixing ratios
    # is an arbitrary function.
    theta = 0.5 * np.pi
    ang_cor = AngularCorrelation(
        State(0, POSITIVE),
        [
            [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.5), State(2, NEGATIVE)],
            [Transition(ELECTRIC, 2, MAGNETIC, 4, -0.5), State(0, POSITIVE)],
        ],
    )
    ana_pow = AnalyzingPower(ang_cor)

    assert ana_pow.evaluate(0.5, ["delta", lambda x: -x], theta) == ana_pow(theta)

    # Test AnalyzingPower.evaluate when the input is a numpy array
    ang_cor_matrix_manual = [
        [
            AnalyzingPower(
                AngularCorrelation(
                    State(0, POSITIVE),
                    [
                        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
                        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.1), State(4, POSITIVE)],
                    ],
                )
            )(theta),
            AnalyzingPower(
                AngularCorrelation(
                    State(0, POSITIVE),
                    [
                        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
                        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.2), State(4, POSITIVE)],
                    ],
                )
            )(theta),
        ],
        [
            AnalyzingPower(
                AngularCorrelation(
                    State(0, POSITIVE),
                    [
                        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
                        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.3), State(4, POSITIVE)],
                    ],
                )
            )(theta),
            AnalyzingPower(
                AngularCorrelation(
                    State(0, POSITIVE),
                    [
                        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
                        [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.4), State(4, POSITIVE)],
                    ],
                )
            )(theta),
        ],
    ]

    ang_cor_matrix = AnalyzingPower(
        AngularCorrelation(
            State(0, POSITIVE),
            [
                [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0), State(2, NEGATIVE)],
                [Transition(ELECTRIC, 2, MAGNETIC, 4, 0.1), State(4, POSITIVE)],
            ],
        )
    ).evaluate(np.array([[0.1, 0.2], [0.3, 0.4]]), [0.0, "delta"], theta=theta)

    assert np.allclose(ang_cor_matrix, ang_cor_matrix_manual)
