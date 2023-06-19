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

#include <cmath>
#include <string>
#include <vector>

using std::to_string;
using std::vector;

#include <gsl/gsl_sf.h>

#include "alpaca/UvCoefficient.hh"

namespace alpaca {

double UvCoefficient::phase_norm_6j_symbol(const int two_nu, const int two_j,
                                           const int two_L, const int two_jp) {

  // Definition of Fagg and Hanna \cite FaggHanna1959 [Eq. (I-1') and the
  // expression below that one]. Causes some tests to fail.
  //
  // const int phase_factor = (((two_jp - two_j - two_L)/2) % 2) == 0 ? 1 : -1;

  // return phase_factor
  // *sqrt(
  //     (two_jp + 1)*(two_j + 1)
  // )
  // *gsl_sf_coupling_6j(
  //     two_j, two_j, two_nu,
  //     two_jp, two_jp, two_L);

  // Definition of Biedenharn \cite AjzenbergSelove1960 (Sec. 1.a.1.iii)
  const int phase_factor = (((two_j + two_jp + two_L) / 2) % 2) == 0 ? 1 : -1;

  return phase_factor * sqrt((two_jp + 1) * (two_j + 1)) *
         gsl_sf_coupling_6j(two_j, two_nu, two_j, two_jp, two_L, two_jp);
}

std::string
UvCoefficient::string_representation(const int n_digits,
                                     vector<string> variable_names) const {

  const string delta_variable =
      variable_names.empty() ? "\\delta" : variable_names[0];

  if (n_digits != 0) {
    return float_string_representation(n_digits, value_L) + "+" +
           float_string_representation(n_digits, value_Lp) + "\\times" +
           delta_variable + "^{2}";
  }

  return "U_{" + to_string(two_nu / 2) + "}\\left(" + to_string(two_j / 2) +
         "," + to_string(two_L / 2) + "," + to_string(two_jp / 2) + "\\right)" +
         "+" + "U_{" + to_string(two_nu / 2) + "}\\left(" +
         to_string(two_j / 2) + "," + to_string(two_Lp / 2) + "," +
         to_string(two_jp / 2) + "\\right)" + delta_variable + "^{2}";
}

} // namespace alpaca
