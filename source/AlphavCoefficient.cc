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

#include <sstream>

using std::stringstream;

#include "alpaca/AlphavCoefficient.hh"

namespace alpaca {

string AlphavCoefficient::string_representation(
    const int n_digits, const vector<string> variable_names) const {

  const string multipole_mixing_ratio_variable =
      variable_names.empty() ? "\\delta" : variable_names[0];

  const string times = (n_digits != 0 ? "\\times" : "");

  return string("(-1)") + (n_digits != 0 ? "\\times" : "") +
         constant_kappa_coefficient.string_representation(n_digits, {}) +
         times + constant_f_coefficient.string_representation(n_digits, {}) +
         "+2" + times +
         linear_kappa_coefficient.string_representation(n_digits, {}) + times +
         linear_f_coefficient.string_representation(n_digits, {}) + times +
         multipole_mixing_ratio_variable + "+" +
         quadratic_kappa_coefficient.string_representation(n_digits, {}) +
         times + quadratic_f_coefficient.string_representation(n_digits, {}) +
         times + multipole_mixing_ratio_variable + "^{2}";
}

} // namespace alpaca
