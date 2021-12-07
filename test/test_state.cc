/*
    This file is part of alpaca.

    alpaca is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    alpaca is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with alpaca.  If not, see <https://www.gnu.org/licenses/>.

    Copyright (C) 2021 Udo Friman-Gayer
*/

#include <cassert>

#include "State.hh"

int main(){
    // Test IO of the State class.

    bool error_thrown = false;
    const State state(2);

    // Error: Unknown parities cannot be converted to string.
    // Usually, a user should never have to call this function.
    try{
        state.parity_str_rep(parity_unknown);
    } catch (const std::runtime_error &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Negative angular momentum quantum number.
    try{
        State state(-1);
    } catch (const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Unfortunately, the State class still allows this at the moment:
    // A user may falsely assume that a half-integer spin is given by a floating-point number.
    // In this case, alpaca will quietly convert the false argument to an integer, and the
    // user will end up with a spin of 1 instead of the desired 3/2.

    State state_with_implicitly_converted_spin(1.5);

    // Test alternative constructor that takes an excitation energy as the second argument.
    
    State state_initialized_with_spin_and_energy(1, 1.);

    // Error: Negative excitation energy given.
    try{
        State state(1, -1.);
    } catch (const std::invalid_argument &e) {
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;
}