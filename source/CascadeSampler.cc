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

#include <iostream>

using std::cout;
using std::endl;

#include <memory>

using std::make_shared;

#include "CascadeSampler.hh"
#include "DirectionSampler.hh"

CascadeSampler::CascadeSampler(vector<AngularCorrelation> &cascade,
                               const int seed, const unsigned int max_tri)
    : initial_direction_random(true), PhiThetaPsi({0., 0., 0.}),
      uniform_direction_sampler(
          []([[maybe_unused]] const double theta,
             [[maybe_unused]] const double phi) { return 1.; },
          1., seed, max_tri) {
  for (size_t i = 0; i < cascade.size(); ++i) {
    angular_correlation_samplers.push_back(
        make_shared<AngCorrRejectionSampler>(cascade[i], seed, max_tri));
  }
}

CascadeSampler::CascadeSampler(vector<AngularCorrelation> &cascade,
                               const int seed,
                               const array<double, 3> PhiThetaPsi,
                               const unsigned int max_tri)
    : initial_direction_random(false), PhiThetaPsi(PhiThetaPsi),
      uniform_direction_sampler(
          []([[maybe_unused]] const double theta,
             [[maybe_unused]] const double phi) { return 1.; },
          1., seed, max_tri){
  for (size_t i = 0; i < cascade.size(); ++i) {
    angular_correlation_samplers.push_back(
        make_shared<AngCorrRejectionSampler>(cascade[i], seed, max_tri));
  }
}

CascadeSampler::CascadeSampler(vector<shared_ptr<DirectionSampler>> cascade,
                               const int seed,
                               const array<double, 3> PhiThetaPsi,
                               const unsigned int max_tri)
    : initial_direction_random(false), PhiThetaPsi(PhiThetaPsi),
      angular_correlation_samplers(cascade),
      uniform_direction_sampler(
          []([[maybe_unused]] const double theta,
             [[maybe_unused]] const double phi) { return 1.; },
          1., seed, max_tri){}

vector<array<double, 2>> CascadeSampler::operator()() {
  vector<array<double, 2>> directions;
  vector<array<double, 3>> reference_frames;

  if (initial_direction_random) {
    const array<double, 2> initial_theta_phi = uniform_direction_sampler();
    reference_frames.push_back(
        {0., initial_theta_phi[0], phi_to_Psi(initial_theta_phi[1])});
  } else {
    reference_frames.push_back(PhiThetaPsi);
  }

  directions.push_back(
      {reference_frames[0][1], Psi_to_phi(reference_frames[0][2])});

  array<double, 2> theta_phi_random;
  for (size_t i = 0; i < angular_correlation_samplers.size(); ++i) {
    theta_phi_random =
        angular_correlation_samplers[i]->operator()(reference_frames[i]);
    directions.push_back(theta_phi_random);
    reference_frames.push_back({0., theta_phi_random[0], 0.});
  }

  return directions;
}
