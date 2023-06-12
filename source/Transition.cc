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

#include <string>

#include "alpaca/Transition.hh"

using std::string;
using std::to_string;

namespace alpaca {

string Transition::str_rep(const State initial_state,
                           const State final_state) const {

  string string_representation = initial_state.str_rep() + " -- ( ";

  if (em_char != EMCharacter::unknown) {
    string_representation += em_str_rep(em_char) + to_string(two_L / 2);
  } else {
    string_representation += to_string(two_L / 2);
  }

  string_representation += " , ";

  if (em_charp != EMCharacter::unknown) {
    string_representation += em_str_rep(em_charp) + to_string(two_Lp / 2);
  } else {
    string_representation += to_string(two_Lp / 2);
  }

  string_representation += " ) --> " + final_state.str_rep();

  return string_representation;
}

} // namespace alpaca
