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

#include <array>

using std::array;

#include <cassert>

#include <gsl/gsl_math.h>

#include "EulerAngleRotation.hh"
#include "SphereRejectionSampler.hh"
#include "TestUtilities.hh"

/**
 * This script tests the correctness of SphereRejectionSampler by sampling from
 * a distribution that is 1 for \f$\varphi < \pi\f$ and 0 otherwise. For an
 * optimum choice \f$W_\mathrm{max} = 1\f$, this should result in an efficiency
 * of \f$\epsilon = 0.5\f$. For the choice \f$W_\mathrm{max} = 2\f$, this should
 * result in an efficiency of \f$\epsilon = 0.25\f$.
 */
int main() {

  SphereRejectionSampler sph_rej_sam(
      []([[maybe_unused]] const double theta, const double phi) {
        return phi < M_PI ? 1. : 0.;
      },
      1., 0);

  double efficiency = sph_rej_sam.estimate_efficiency(1e5);

  test_numerical_equality<double>(efficiency, 0.5, 1e-3);

  SphereRejectionSampler sph_rej_sam_2(
      []([[maybe_unused]] const double theta, const double phi) {
        return phi < M_PI ? 1. : 0.;
      },
      2., 0);

  efficiency = sph_rej_sam_2.estimate_efficiency(1e5);

  test_numerical_equality<double>(efficiency, 0.25, 5e-3);

  // Test the case in which no vector can be found.
  SphereRejectionSampler sph_rej_sam_3(
      []([[maybe_unused]] const double theta,
         [[maybe_unused]] const double phi) { return -1.; },
      0.5, 0);
  const pair<unsigned int, array<double, 3>> theta_phi_default =
      sph_rej_sam_3.sample();
  assert(theta_phi_default.first == 1000);
  assert(theta_phi_default.second[0] == 0.);
  assert(theta_phi_default.second[1] == 0.);
  assert(theta_phi_default.second[2] == 0.);
}