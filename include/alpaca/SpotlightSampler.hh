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

#pragma once

#include <random>

#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/ReferenceFrameSampler.hh"

using std::mt19937;
using std::uniform_real_distribution;

namespace alpaca {

class SpotlightSampler : public ReferenceFrameSampler {
public:
  SpotlightSampler(const CoordDir a_theta_phi, const double a_opening_angle,
                   const unsigned long a_seed)
      : theta_phi(a_theta_phi), opening_angle(a_opening_angle),
        u_min(0.5 * (1. + cos(opening_angle))), seed(a_seed),
        random_engine(mt19937(seed)) {}
  SpotlightSampler(const CoordDir a_theta_phi, const unsigned long a_seed)
      : SpotlightSampler(a_theta_phi, 0., a_seed) {}
  SpotlightSampler(const CoordDir a_theta_phi, const double distance,
                   const double radius, const unsigned long a_seed)
      : SpotlightSampler(a_theta_phi, asin(radius / distance), a_seed) {}

  pair<unsigned int, EulerAngles> sample() override;

protected:
  CoordDir theta_phi;
  double opening_angle;
  double u_min{0.5};
  unsigned long seed;

  mt19937 random_engine; /**< Deterministic random number engine. */
  uniform_real_distribution<double>
      uniform_random; /**< Uniform distribution from which all random numbers
                         are derived here. */
};

} // namespace alpaca
