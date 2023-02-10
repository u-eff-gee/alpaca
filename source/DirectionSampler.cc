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

#include <utility>

using std::pair;

#include "DirectionSampler.hh"

array<double, 2> DirectionSampler::operator()() {
  pair<unsigned int, array<double, 2>> sampled_theta_phi = sample();

  return {sampled_theta_phi.second[0], sampled_theta_phi.second[1]};
}

array<double, 2>
DirectionSampler::operator()(const array<double, 3> euler_angles) {
  return euler_angle_rotation.rotate(operator()(), euler_angles);
}