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
#include <gsl/gsl_math.h>

#include "alpaca/AngularCorrelation.hh"
#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/State.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/Transition.hh"
#include "alpaca/W_dir_dir.hh"
#include "alpaca/W_pol_dir.hh"

using alpaca::AngularCorrelation;
using alpaca::EMCharacter;
using alpaca::Parity;
using alpaca::State;
using alpaca::test_numerical_equality;
using alpaca::Transition;
using alpaca::W_dir_dir;
using alpaca::W_pol_dir;

/**
 * Test whether AngularCorrelation maps to the correct correlation functions and
 * whether the rotation works.
 *
 * The rotation is tested by comparing correlations with polarized incident
 * gamma rays for different electromagnetic (EM) characters of the first
 * transition. The difference in the EM character leads to a rotation of the
 * distribution by \f$\pi / 2\f$ around the z axis.
 */
int main() {

  const double epsilon = 1e-8;

  // Test unpolarized angular correlation
  AngularCorrelation ang_corr_0_1_0{
      State(0, Parity::unknown),
      {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
        State(2, Parity::negative)},
       {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
        State(0, Parity::unknown)}}};

  W_dir_dir w_dir_dir_0_1_0{
      State(0, Parity::positive),
      {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
        State(2, Parity::negative)},
       {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
        State(0, Parity::positive)}}};

  for (double theta = 0.; theta < M_PI; theta += 0.5) {
    for (double phi = 0.; phi < M_2_PI; phi += 0.5) {

      test_numerical_equality<double>(ang_corr_0_1_0(theta, phi),
                                      w_dir_dir_0_1_0(theta, phi), epsilon);
    }
  }

  // Test polarized angular correlation and rotation
  AngularCorrelation ang_corr_0p_1p_0p{
      State(0, Parity::positive),
      {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
        State(2, Parity::positive)},
       {Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
        State(0, Parity::positive)}}};

  W_pol_dir w_pol_dir_0p_1p_0p{
      State(0, Parity::positive),
      {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
        State(2, Parity::positive)},
       {Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
        State(0, Parity::positive)}}};

  W_pol_dir w_pol_dir_0p_1m_0p{
      State(0, Parity::positive),
      {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
        State(2, Parity::negative)},
       {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
        State(0, Parity::positive)}}};

  for (double theta = 0.; theta < M_PI; theta += 0.5) {
    for (double phi = 0.; phi < M_2_PI; phi += 0.5) {

      test_numerical_equality<double>(ang_corr_0p_1p_0p(theta, phi),
                                      w_pol_dir_0p_1p_0p(theta, phi), epsilon);
      // This double-loop test was originally intended to test the capability of
      // angular correlations to rotate themselves.
      // Since explicit rotations are preferred now, this at least tests the
      // symmetry of the implemented distributions. At the moment, rotation is
      // achieved using a trick.
      test_numerical_equality<double>(ang_corr_0p_1p_0p(theta, phi),
                                      w_pol_dir_0p_1m_0p(theta, phi - M_PI_2),
                                      epsilon);
    }
  }

  // Test the copy constructor
  AngularCorrelation ang_corr_0_1_0_prime = ang_corr_0_1_0;
  assert(ang_corr_0_1_0_prime(0.1, 0.2) == ang_corr_0_1_0(0.1, 0.2));
}
