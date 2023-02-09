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

#include "SpherePointSampler.hh"
#include "TestUtilities.hh"

int main() {
  // Test the method that samples points in cartesian coordinates by verifying
  // that they lie on the surface of a sphere.
  SpherePointSampler sph_poi_sam;
  const double radius = 2.;
  const double radius_squared = radius * radius;
  const unsigned int n = 100;
  array<vector<double>, 3> x_y_z = sph_poi_sam.sample_cartesian(n, radius);
  for (size_t i = 0; i < n; ++i) {
    test_numerical_equality<double>(x_y_z[0][i] * x_y_z[0][i] +
                                        x_y_z[1][i] * x_y_z[1][i] +
                                        x_y_z[2][i] * x_y_z[2][i],
                                    radius_squared, 1e-5);
  }

  // Test various error messages that should be displayed when the numerical
  // fixed-point searches fail.
  bool error_thrown = false;

  try {
    sph_poi_sam.find_c(2, 1e-8, 2);
  } catch (const runtime_error &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  try {
    sph_poi_sam.find_Theta_j(1, 1, 0.1, 1e-8, 2);
  } catch (const runtime_error &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Test the elliptic integral of the 1st kind with arbitrary m.
  // After testing the alpaca code for months, it turned out that the code in
  // which the elliptic integrals of the first kind are called for arbitrary phi
  // and m is actually never used. The test below compares to a value that was
  // calculated using Mathematica's \cite Weisstein2021 EllipticF(0.1, -1)
  test_numerical_equality<double>(
      sph_poi_sam.elliptic_integral_1st_kind_arbitrary_m(0.1, -0.1), 0.099,
      1e-3);
}