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
#include <cmath>

#include <gsl/gsl_math.h>

#include "alpaca/AvCoefficient.hh"
#include "alpaca/State.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/Transition.hh"
#include "alpaca/UvCoefficient.hh"
#include "alpaca/W_dir_dir.hh"

using alpaca::AvCoefficient;
using alpaca::State;
using alpaca::test_numerical_equality;
using alpaca::Transition;
using alpaca::UvCoefficient;
using alpaca::W_dir_dir;

/**
 * \brief Analytical expression for the angular distribution of a \f$0 \to 2 \to
 * 0\f$ sequence
 *
 * The expression is { Eq. (11) in \cite Kneissl1996}:
 *
 * \f[
 * 		W\left( \theta \right) = \frac{5}{4} \left[ 1 - 3\cos^2 \left(
 * \theta \right) + 4\cos^4 \left( \theta \right) \right]. \f]
 */
double w_dir_dir_0_2_0(const double theta) {
  return (1.25 - 3.75 * pow(cos(theta), 2) + 5. * pow(cos(theta), 4));
}

int main() {

  const double epsilon{1e-7};

  double w_dir_dir_num{0.};
  double w_dir_dir_ana{0.};

  W_dir_dir w_dir_dir_2(State(0), {
                                      {Transition(4, 6, 0.), State(4)},
                                      {Transition(4, 6, 0.), State(0)},
                                  });

  for (double theta = 0.; theta < M_PI; theta += 0.5) {
    w_dir_dir_num = w_dir_dir_2(theta);
    w_dir_dir_ana = w_dir_dir_0_2_0(theta);

    test_numerical_equality<double>(w_dir_dir_num, w_dir_dir_ana, epsilon);
  }

  // Test string representation.
  // As a test case, use the 0->1->2 direction-direction correlation in
  // Sec. "4 Numerical example" of Ref. \cite Iliadis2021.
  W_dir_dir w_0_1_2(State(0), {
                                  {Transition(2, 4, 0.), State(2)},
                                  {Transition(2, 4, 0.), State(4)},
                              });

  const string str_rep =
      string("\\left[") +
      AvCoefficient(0, 2, 4, 0, 2).string_representation(0, {"\\delta_1"}) +
      "\\right]\\\\" + "\\times\\left[" +
      AvCoefficient(0, 2, 4, 4, 2).string_representation(0, {"\\delta_2"}) +
      "\\right]\\\\" +
      "\\times P_{0}\\left[\\cos\\left(\\theta\\right)\\right]\\\\" +
      "+\\left[" +
      AvCoefficient(4, 2, 4, 0, 2).string_representation(0, {"\\delta_1"}) +
      "\\right]\\\\" + "\\times\\left[" +
      AvCoefficient(4, 2, 4, 4, 2).string_representation(0, {"\\delta_2"}) +
      "\\right]\\\\" +
      "\\times P_{2}\\left[\\cos\\left(\\theta\\right)\\right]";
  assert(w_0_1_2.string_representation() == str_rep);

  // Test string representation for cascade with intermediate step.
  W_dir_dir w_0_1_1_2(State(0), {
                                    {Transition(2, 4, 0.), State(2)},
                                    {Transition(2, 4, 0.), State(2)},
                                    {Transition(2, 4, 0.), State(4)},
                                });

  const string str_rep_2 =
      string("\\left[") +
      AvCoefficient(0, 2, 4, 0, 2).string_representation(0, {"\\delta_1"}) +
      "\\right]\\\\" + "\\times\\left[" +
      UvCoefficient(0, 2, 2, 4, 0., 2).string_representation(0, {"\\delta_2"}) +
      "\\right]\\\\" + "\\times\\left[" +
      AvCoefficient(0, 2, 4, 4, 2).string_representation(0, {"\\delta_3"}) +
      "\\right]\\\\" +
      "\\times P_{0}\\left[\\cos\\left(\\theta\\right)\\right]\\\\" +
      "+\\left[" +
      AvCoefficient(4, 2, 4, 0, 2).string_representation(0, {"\\delta_1"}) +
      "\\right]\\\\" + "\\times\\left[" +
      UvCoefficient(4, 2, 2, 4, 0., 2).string_representation(0, {"\\delta_2"}) +
      "\\right]\\\\" + "\\times\\left[" +
      AvCoefficient(4, 2, 4, 4, 2).string_representation(0, {"\\delta_3"}) +
      "\\right]\\\\" +
      "\\times P_{2}\\left[\\cos\\left(\\theta\\right)\\right]";
  assert(w_0_1_1_2.string_representation() == str_rep_2);
}
