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

#include "AvCoefficient.hh"

AvCoefficient::AvCoefficient(const int two_nu, const int two_L,
                             const int two_Lp, const int two_jn,
                             const int two_j)
    : two_nu(two_nu), two_L(two_L), two_Lp(two_Lp), two_jn(two_jn),
      two_j(two_j), constant_f_coefficient(two_nu, two_L, two_L, two_jn, two_j),
      linear_f_coefficient(two_nu, two_L, two_Lp, two_jn, two_j),
      quadratic_f_coefficient(two_nu, two_Lp, two_Lp, two_jn, two_j) {
  constant_coefficient = constant_f_coefficient.get_value();
  linear_coefficient = 2. * linear_f_coefficient.get_value();
  quadratic_coefficient = quadratic_f_coefficient.get_value();
}

double AvCoefficient::operator()(const double delta) const {

  return constant_coefficient + delta * linear_coefficient +
         delta * delta * quadratic_coefficient;
}

string AvCoefficient::string_representation(
    const unsigned int n_digits, const vector<string> variable_names) const {

  string multipole_mixing_ratio_variable =
      variable_names.size() ? variable_names[0] : "\\delta";

  return constant_f_coefficient.string_representation(n_digits, {}) + "+" +
         "2" + (n_digits ? "\\times" : "") +
         linear_f_coefficient.string_representation(n_digits, {}) +
         (n_digits ? "\\times" : "") + multipole_mixing_ratio_variable + "+" +
         quadratic_f_coefficient.string_representation(n_digits, {}) +
         (n_digits ? "\\times" : "") + multipole_mixing_ratio_variable + "^{2}";
}