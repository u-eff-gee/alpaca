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

#include "alpaca/AlphavCoefficient.hh"
#include "alpaca/EvCoefficient.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/Transition.hh"

using alpaca::AlphavCoefficient;
using alpaca::EMCharacter;
using alpaca::EvCoefficient;
using alpaca::test_numerical_equality;

int main() {

  // Test the E_v coefficients of Biedenharn \cite AjzenbergSelove1960 against
  // the alpha_v coefficients of Fagg and Hanna \cite FaggHanna1959. Both
  // references introduce a sign into the expression for the
  // polarization-direction correlation which is determined by the
  // electromagnetic character of the associated radiation. Note that in Ref.
  // \cite AjzenbergSelove1960 [Eq. (85a)], this sign is part of the E_v
  // coefficients, while it is factored out in Ref. \cite FaggHanna1959 [Eq.
  // (I-8)].
  AlphavCoefficient alphav_coef(4, 4, 6, 7, 3);
  EvCoefficient ev_coef(4, EMCharacter::electric, 4, EMCharacter::magnetic, 6,
                        7, 3);
  const double epsilon = 1e-4;

  double ev_coef_num{0.};

  // Eq. (88b) in \cite AjzenbergSelove1960
  ev_coef_num = ev_coef(0.);
  test_numerical_equality<double>(ev_coef_num, 0.07143, epsilon);

  // The primary transition is an electric transition.
  // Equation (85a) in \cite AjzenbergSelove1960 shows that the sign of the term
  // which is independent of delta is Parity::positive. Equations (I-8) and
  // (I-9) in \cite FaggHanna1959 show that the sign of the same term is
  // Parity::positive as well, i.e. no additional sign is needed here.
  //
  // Note that the notation for the sign in Ref. \cite FaggHanna1959 may easily
  // be misunderstood, because it is not stated clearly. However, a cross check
  // with other literature, for example Ref. \cite Kneissl1996, shows that
  // indeed, 'the + or - sign is to be employed according as the ... is
  // electric or magnetic in character' \cite FaggHanna1959 means
  // that 'electric' and + belong together.
  test_numerical_equality<double>(ev_coef_num, alphav_coef(0.), epsilon);

  // Eq. (89b) in \cite AjzenbergSelove1960
  EvCoefficient ev_coef_2(4, EMCharacter::magnetic, 2, EMCharacter::electric, 4,
                          3, 3);
  ev_coef_num = ev_coef_2(1.);
  test_numerical_equality<double>(ev_coef_num, 0.20000 - 0.25820, epsilon);

  AlphavCoefficient alphav_coef_2(4, 2, 4, 3, 3);
  // Here, the additional sign is needed.
  test_numerical_equality<double>(ev_coef_num, -alphav_coef_2(1.), epsilon);
}
