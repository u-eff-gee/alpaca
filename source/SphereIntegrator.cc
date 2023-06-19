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

#include <vector>

#include <gsl/gsl_math.h>

#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/SphereIntegrator.hh"

using std::vector;

namespace alpaca {

double SphereIntegrator::operator()(double f(double theta, double phi),
                                    const unsigned int n,
                                    bool is_in_omega(double theta,
                                                     double phi)) {

  vector<CoordDir> theta_phi = sph_poi_samp.sample(n);

  double integral = 0.;

  for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
    if (is_in_omega(theta_phi[i][0], theta_phi[i][1])) {
      integral += f(theta_phi[i][0], theta_phi[i][1]);
    }
  }

  return 4. * M_PI / static_cast<double>(n) * integral;
}

} // namespace alpaca
