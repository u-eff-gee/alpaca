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
#include "EulerAngleRotation.hh"
#include "ReferenceFrameSampler.hh"

CascadeSampler::CascadeSampler(
    vector<shared_ptr<ReferenceFrameSampler>> cascade)
    : angular_correlation_samplers(cascade) {}

vector<array<double, 3>> CascadeSampler::operator()() {
  vector<array<double, 3>> reference_frames(
      angular_correlation_samplers.size());

  reference_frames[0] = angular_correlation_samplers[0]->operator()();

  for (size_t i = 1; i < angular_correlation_samplers.size(); ++i) {
    reference_frames[i] = euler_angle_transform::rotate(
        angular_correlation_samplers[i]->operator()(), reference_frames[i - 1]);
  }

  return reference_frames;
}
