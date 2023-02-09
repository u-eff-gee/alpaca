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

#include "CascadeRejectionSampler.hh"

CascadeRejectionSampler::CascadeRejectionSampler(
    vector<AngularCorrelation> &cascade, const int seed,
    const unsigned int max_tri)
    : initial_direction_random(true), PhiThetaPsi({0., 0., 0.}),
      return_first_direction(true),
      uniform_direction_sampler(
          []([[maybe_unused]] const double theta,
             [[maybe_unused]] const double phi) { return 1.; },
          1., seed, max_tri),
      euler_angle_rotation(EulerAngleRotation()) {
  for (size_t i = 0; i < cascade.size(); ++i) {
    angular_correlation_samplers.push_back(
        AngCorrRejectionSampler(cascade[i], seed, max_tri));
  }
}

CascadeRejectionSampler::CascadeRejectionSampler(
    vector<AngularCorrelation> &cascade, const int seed,
    const array<double, 3> PhiThetaPsi, const bool return_first_direction,
    const unsigned int max_tri)
    : initial_direction_random(false), PhiThetaPsi(PhiThetaPsi),
      return_first_direction(return_first_direction),
      uniform_direction_sampler(
          []([[maybe_unused]] const double theta,
             [[maybe_unused]] const double phi) { return 1.; },
          1., seed, max_tri),
      euler_angle_rotation(EulerAngleRotation()) {
  for (size_t i = 0; i < cascade.size(); ++i) {
    angular_correlation_samplers.push_back(
        AngCorrRejectionSampler(cascade[i], seed, max_tri));
  }
}

vector<array<double, 2>> CascadeRejectionSampler::operator()() {
  vector<array<double, 2>> directions;
  vector<array<double, 3>> reference_frames;

  if (initial_direction_random) {
    const array<double, 2> initial_theta_phi = uniform_direction_sampler();
    reference_frames.push_back(
        {0., initial_theta_phi[0], phi_to_Psi(initial_theta_phi[1])});
  } else {
    reference_frames.push_back(PhiThetaPsi);
  }

  if (return_first_direction) {
    directions.push_back(
        {reference_frames[0][1], Psi_to_phi(reference_frames[0][2])});
  }

  array<double, 2> theta_phi_random;
  for (size_t i = 0; i < angular_correlation_samplers.size(); ++i) {
    theta_phi_random = angular_correlation_samplers[i](reference_frames[i]);
    directions.push_back(theta_phi_random);
    reference_frames.push_back(
        {0., theta_phi_random[0], phi_to_Psi(theta_phi_random[1])});
  }

  return directions;
}
