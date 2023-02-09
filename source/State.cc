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

    Copyright (C) 2021-2023 Udo Friman-Gayer
*/

#include "State.hh"

string State::parity_str_rep(const Parity parity) const {
		
    if(parity == positive){
        return "+";
    } 
    if(parity == negative){
        return "-";
    }

    throw runtime_error("No string representation for unknown parity.");

}

string State::spin_str_rep(const int two_J) const {

    if(two_J % 2 == 0){
        return to_string(two_J/2);
    }

    return to_string(two_J) + "/2";
    
}

string State::str_rep() const {

    if(parity != parity_unknown){
        return spin_str_rep(two_J) + "^" + parity_str_rep(parity); 
    }

    return spin_str_rep(two_J);
}

int State::check_two_J(const int two_J) const {
    
    if(two_J < 0){
        throw invalid_argument("two_J must be a nonnegative integer.");
    }

    return two_J;
}

double State::check_excitation_energy(const double e_x) const {
    
    if(e_x < 0.){
        throw invalid_argument("Excitation energy must not be negative.");
    }

    return e_x;
}