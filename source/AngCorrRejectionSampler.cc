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

#include "AngCorrRejectionSampler.hh"
#include "EulerAngleRotation.hh"

AngCorrRejectionSampler::AngCorrRejectionSampler(AngularCorrelation &w,
                                                 const int seed,
                                                 const unsigned int max_tri)
    : SphereRejectionSampler(nullptr, w.get_upper_limit(), seed, max_tri),
      angular_correlation(w.get_initial_state(), w.get_cascade_steps()) {}

pair<unsigned int, array<double, 3>> AngCorrRejectionSampler::sample() {

  array<double, 2> theta_phi;
  double dis_val;

  for (unsigned int i = 0; i < max_tries; ++i) {

    theta_phi = sample_theta_phi();
    dis_val = uniform_random(random_engine) * distribution_maximum;

    if (dis_val <= angular_correlation(theta_phi[0], theta_phi[1])) {
      return {i + 1, euler_angle_transform::from_spherical(
                         theta_phi, uniform_random(random_engine) * 2. * M_PI)};
    }
  }

  return {max_tries, {0., 0., 0.}};
}