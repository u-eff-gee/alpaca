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

#include <gsl/gsl_sf.h>

#include "alpaca/State.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/Transition.hh"
#include "alpaca/W_dir_dir.hh"

using alpaca::EMCharacter;
using alpaca::Parity;
using alpaca::State;
using alpaca::test_numerical_equality;
using alpaca::Transition;
using alpaca::W_dir_dir;

/**
 * Eq. (68) in Ref. \cite AjzenbergSelove1960.
 */
double w_dir_dir_6_4_3_1(const double theta) {
  return 1. + 0.10204 * gsl_sf_legendre_Pl(2, cos(theta)) +
         0.00907 * gsl_sf_legendre_Pl(4, cos(theta));
}

/**
 * \brief Test angular correlations with unobserved intermediate transitions.
 *
 * The present test uses an example given in Sec. 1.a.1.iii [Eqs. (66) - (68)]
 * of Ref. \cite AjzenbergSelove1960.
 */
int main() {

  const double epsilon = 1e-4;

  double w_num{0.}, w_ana{0.};

  W_dir_dir w_dir_dir(
      State(12, Parity::unknown),
      {
          {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
           State(8, Parity::unknown)},
          {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
           State(6, Parity::unknown)},
          {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
           State(2, Parity::unknown)},
      });

  for (double theta = 0.; theta < M_PI; theta += 0.5) {
    for (double phi = 0.; phi < M_2_PI; phi += 0.5) {

      w_num = w_dir_dir(theta);
      w_ana = w_dir_dir_6_4_3_1(theta);

      test_numerical_equality<double>(w_num, w_ana, epsilon);
    }
  }
}
