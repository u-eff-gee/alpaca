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

#include <cassert>
#include <memory>

using std::make_unique;
using std::unique_ptr;

#include "alpaca/TestUtilities.hh"
#include "alpaca/UvCoefficient.hh"

using alpaca::test_numerical_equality;
using alpaca::UvCoefficient;

int main() {

  // Test the calculation of \f$U_\nu\f$ coefficients using the example given in
  // Sec. 1.a.1.iii [Eqs. (66) - (68)] of Ref. \cite AjzenbergSelove1960, which
  // gives the values of two coefficients explicitly.

  unique_ptr<UvCoefficient> uv_coef = make_unique<UvCoefficient>(0, 0, 0, 0);
  const double epsilon = 1e-4;

  // Test the coefficient \f$U_2\f$
  uv_coef = make_unique<UvCoefficient>(4, 8, 2, 6);
  test_numerical_equality<double>(uv_coef->get_value(), 0.90469, epsilon);

  // Test the coefficient \f$U_4\f$
  uv_coef = make_unique<UvCoefficient>(8, 8, 2, 6);
  test_numerical_equality<double>(uv_coef->get_value(), 0.68138, epsilon);

  // Test the string representation for both constructors.
  assert(uv_coef->string_representation() ==
         "U_{4}\\left(4,1,3\\right)+U_{4}\\left(4,2,3\\right)\\delta^{2}");
  assert(uv_coef->string_representation(3) == "0.681+0\\times\\delta^{2}");

  uv_coef = make_unique<UvCoefficient>(8, 8, 2, 4, 0., 6);
  assert(uv_coef->string_representation() ==
         "U_{4}\\left(4,1,3\\right)+U_{4}\\left(4,2,3\\right)\\delta^{2}");
  assert(uv_coef->string_representation(3) == "0.681+0\\times\\delta^{2}");
}
