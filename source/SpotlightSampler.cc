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

#include <cmath>
#include <utility>

#include <gsl/gsl_math.h>

#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/SpotlightSampler.hh"

using std::pair;

using alpaca::CoordDir;
using alpaca::EulerAngles;

namespace alpaca {

pair<unsigned int, EulerAngles> SpotlightSampler::sample() {
  if (opening_angle == 0.0) {
    return {1, euler_angle_transform::from_spherical(theta_phi)};
  }

  const double theta =
      acos(2.0 * (u_min + (1.0 - u_min) * uniform_random(random_engine)) - 1.0);
  const double phi = 2.0 * M_PI * uniform_random(random_engine);

  return {1, euler_angle_transform::from_spherical({theta, phi}, 0.)};
}

} // namespace alpaca
