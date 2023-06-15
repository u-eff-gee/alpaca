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
#    Copyright (C) 2021-2023 Udo Friman-Gayer

import pytest

from alpaca import EMCharacter, Transition


def test_transition():
    # Multipolarities must be given with a factor of two, i.e. there can be no odd values.
    with pytest.raises(ValueError):
        Transition(EMCharacter.electric, 1, EMCharacter.magnetic, 4, 0.0)
    with pytest.raises(ValueError):
        Transition(EMCharacter.electric, 2, EMCharacter.magnetic, 3, 0.0)

    # assert Transition(EMCharacter.electric, 2, EMCharacter.magnetic, 4, 0.0).str() == "E1 (M2)"
    # assert Transition(EMCharacter.electric, 2, EMCharacter.magnetic, 4, 0.0).__str__() == "E1 (M2)"
    # assert Transition(EMCharacter.electric, 2, EMCharacter.magnetic, 4, 0.0).str(separator=",") == "E1,(M2)"
    # assert (
    #     Transition(EMCharacter.electric, 2, EMCharacter.magnetic, 4, 0.0).str(secondary_in_parentheses=False)
    #     == "E1 M2"
    # )
    # assert (
    #     Transition(EMCharacter.electric, 2, EMCharacter.magnetic, 4, 0.0).str(always_show_secondary=False)
    #     == "E1"
    # )
    # assert (
    #     Transition(EMCharacter.unknown, 2, EMCharacter.unknown, 4, 0.0).str(always_show_secondary=False)
    #     == "σ1"
    # )
