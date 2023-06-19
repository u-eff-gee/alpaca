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

#pragma once

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using std::string;
using std::stringstream;
using std::vector;

namespace alpaca {

/**'
 * \brief Abstract class for string-representable expressions
 */
class StringRepresentable {

public:
  virtual ~StringRepresentable() = default;
  StringRepresentable() = default;
  StringRepresentable(const StringRepresentable &) = default;
  StringRepresentable &operator=(const StringRepresentable &) = default;
  StringRepresentable(StringRepresentable &&) = default;
  StringRepresentable &operator=(StringRepresentable &&) = default;

  /**
   * \brief Return string representation of expression.
   *
   * This function has an argument n_digits that indicates whether coefficients
   * should be displayed as variables or evaluated numerically. The numerical
   * expressions will be formatted using std::setprecision(n_digits).
   *
   * \param n_digits Determines whether the expression
   * should be evaluated numerically, (n_digits > 0). If yes, indicates how many
   * digits should be displayed (default: 0).
   * \param variable_names Names for the variables of a function (default: {}
   * i.e. use default names).
   *
   * \return String representation.
   */
  [[nodiscard]] virtual string
  string_representation(const int n_digits = 0,
                        const vector<string> variable_names = {}) const = 0;

protected:
  [[nodiscard]] static string float_string_representation(const int n_digits,
                                                          const double number) {
    stringstream str_rep;
    if (number < 0.) {
      str_rep << "\\left(" << std::setprecision(n_digits) << number
              << "\\right)";
    } else {
      str_rep << std::setprecision(n_digits) << number;
    }
    return str_rep.str();
  };
};

} // namespace alpaca
