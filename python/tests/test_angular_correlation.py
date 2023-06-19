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

# As stated in the README: Please note that this test does not ensure that the angular-correlation
# formalism is implemented correctly.
# This is already done in the tests of the C++ code.
# The purpose of this test is to ensure that the python API works correctly.

from alpaca import AngularCorrelation, Parity, State, EMCharacter, Transition


def test_angular_correlation():
    initial_state = State(0, Parity.positive)
    cascade_steps = [
        [
            Transition(EMCharacter.magnetic, 2, EMCharacter.electric, 4, 0.0),
            State(2, Parity.positive),
        ],
        [
            Transition(EMCharacter.magnetic, 2, EMCharacter.electric, 4, 0.5),
            State(4, Parity.positive),
        ],
    ]
    assert len(cascade_steps) == 2
    ang_cor = AngularCorrelation(initial_state, cascade_steps)

    # assert ang_cor(0.1, 0.1) == angular_correlation(
    #     0.1, 0.1, initial_state, cascade_steps
    # )
    # assert ang_cor(0.1, 0.1, Phi_Theta_Psi=(0.1, 0.1, 0.1)) == angular_correlation(
    #     0.1, 0.1, initial_state, cascade_steps, Phi_Theta_Psi=(0.1, 0.1, 0.1)
    # )
