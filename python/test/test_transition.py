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

from alpaca.transition import ELECTRIC, EM_UNKNOWN, MAGNETIC, Transition
import pytest


def test_transition():
    # Multipolarities must be given with a factor of two, i.e. there can be no odd values.
    with pytest.raises(ValueError):
        Transition(ELECTRIC, 1, MAGNETIC, 4, 0.0)
    with pytest.raises(ValueError):
        Transition(ELECTRIC, 2, MAGNETIC, 3, 0.0)

    assert Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0).str() == "E1 (M2)"
    assert Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0).__str__() == "E1 (M2)"
    assert Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0).str(separator=",") == "E1,(M2)"
    assert (
        Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0).str(secondary_in_parentheses=False)
        == "E1 M2"
    )
    assert (
        Transition(ELECTRIC, 2, MAGNETIC, 4, 0.0).str(always_show_secondary=False)
        == "E1"
    )
    assert (
        Transition(EM_UNKNOWN, 2, EM_UNKNOWN, 4, 0.0).str(always_show_secondary=False)
        == "σ1"
    )
