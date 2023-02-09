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

#include "Transition.hh"

Transition::Transition(const int t_L, const int t_Lp, const double del)
    : em_char(em_unknown), two_L(check_two_L(t_L)), em_charp(em_unknown),
      two_Lp(check_two_L(t_Lp)), delta(del) {
  if (two_L == two_Lp) {
    throw invalid_argument(
        "The two multipolarities for a transition may not be equal. This holds "
        "even if the coupling allows only a single multipolarity.");
  }
}

Transition::Transition(const EMCharacter em, const int t_L,
                       const EMCharacter emp, const int t_Lp, const double del)
    : em_char(em), two_L(check_two_L(t_L)), em_charp(emp),
      two_Lp(check_two_L(t_Lp)), delta(del) {
  if (two_L == two_Lp) {
    throw invalid_argument(
        "The two multipolarities for a transition may not be equal. This holds "
        "even if the coupling allows only a single multipolarity.");
  }
}

string Transition::str_rep(const State initial_state,
                           const State final_state) const {

  string string_representation = initial_state.str_rep() + " -- ( ";

  if (em_char != em_unknown) {
    string_representation += em_str_rep(em_char) + to_string(two_L / 2);
  } else {
    string_representation += to_string(two_L / 2);
  }

  string_representation += " , ";

  if (em_charp != em_unknown) {
    string_representation += em_str_rep(em_charp) + to_string(two_Lp / 2);
  } else {
    string_representation += to_string(two_Lp / 2);
  }

  string_representation += " ) --> " + final_state.str_rep();

  return string_representation;
}

int Transition::check_two_L(const int two_L) const {

  if (two_L < 1) {
    throw invalid_argument(
        "two_L (two_Lp) must be a nonzero, nonnegative integer.");
  }

  return two_L;
}