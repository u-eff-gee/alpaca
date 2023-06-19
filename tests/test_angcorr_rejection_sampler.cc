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

#include <gsl/gsl_sf.h>

#include "alpaca/AngCorrRejectionSampler.hh"
#include "alpaca/AngularCorrelation.hh"
#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/SphereRejectionSampler.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/Transition.hh"

namespace euler_angle_transform = alpaca::euler_angle_transform;

using alpaca::AngCorrRejectionSampler;
using alpaca::AngularCorrelation;
using alpaca::EMCharacter;
using alpaca::Parity;
using alpaca::SphereRejectionSampler;
using alpaca::State;
using alpaca::test_numerical_equality;
using alpaca::Transition;
using alpaca::CoordDir;

/**
 * Test the AngCorrRejectionSampler by verifying that the class, using an
 * AngularCorrelation object, does exactly the same as a SphereRejectionSampler,
 * using an equivalent analytical expression of the angular correlation.
 */
int main() {

  const int seed = 0;

  AngularCorrelation ang_cor(
      State(0, Parity::positive),
      {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
        State(2, Parity::positive)},
       {Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
        State(0, Parity::positive)}});

  AngCorrRejectionSampler ang_cor_sam(ang_cor, seed);

  // Analytical expression for the angular correlation above.
  // See, e.g., Eq. (1) in Ref. \cite Pietralla2001.
  SphereRejectionSampler sph_rej_sam(
      [](const double theta, const double phi) {
        return 1. + 0.5 * (gsl_sf_legendre_Pl(2, cos(theta)) +
                           0.5 * cos(2. * phi) *
                               gsl_sf_legendre_Plm(2, 2, cos(theta)));
      },
      ang_cor.get_upper_limit(), seed);

  CoordDir theta_phi_1;
  CoordDir theta_phi_2;

  for (unsigned int n = 0; n < 10; ++n) {
    theta_phi_1 = euler_angle_transform::to_spherical(ang_cor_sam());
    theta_phi_2 = euler_angle_transform::to_spherical(sph_rej_sam());

    test_numerical_equality<double>(2, theta_phi_1.data(), theta_phi_2.data(),
                                    1e-6);
  }

  // Check that the default values
  //
  // theta = 0
  // phi = pi/2
  //
  // corresponding to Euler angles
  //
  // Phi = 0
  // Theta = 0
  // Psi = 0
  //
  // are returned when AngCorrRejectionSampler can not find a valid vector.
  // In order to test it, the max_tri is set to zero.
  // This way, the random sampling is bypassed.
  AngCorrRejectionSampler ang_cor_sam_2(ang_cor, 0, 0);
  theta_phi_1 = euler_angle_transform::to_spherical(ang_cor_sam_2());
  assert(theta_phi_1[0] == 0.);
  assert(theta_phi_1[1] == M_PI_2);
}
