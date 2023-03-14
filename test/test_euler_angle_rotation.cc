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

#include <gsl/gsl_math.h>

#include "EulerAngleRotation.hh"
#include "TestUtilities.hh"

/**
 * \brief Test rotations of 3D vectors using Euler angles.
 *
 * Test by rotating the three canonical Cartesian axes into each other.
 */
int main() {

  const double epsilon = 1e-8;

  array<double, 3> x_axis{1., 0., 0.};
  array<double, 3> y_axis{0., 1., 0.};
  array<double, 3> z_axis{0., 0., 1.};

  array<double, 3> euler_angles{0., 0., 0.};
  array<double, 3> xp_yp_zp{0., 0., 0.};

  // Rotate x axis into y axis
  // Phi   = -pi/2
  // Theta = 0
  // Psi   = 0
  euler_angles = {-M_PI_2, 0., 0.};

  xp_yp_zp = euler_angle_transform::rotate(x_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), y_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, euler_angle_transform::rotate_back(xp_yp_zp, euler_angles).data(),
      x_axis.data(), epsilon);

  // Rotate x axis into z axis
  // Phi   = pi/2
  // Theta = pi/2
  // Psi   = 0
  euler_angles = {M_PI_2, M_PI_2, 0.};

  xp_yp_zp = euler_angle_transform::rotate(x_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), z_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, euler_angle_transform::rotate_back(xp_yp_zp, euler_angles).data(),
      x_axis.data(), epsilon);

  // Rotate y axis into x axis
  // Phi   = pi/2
  // Theta = 0
  // Psi   = 0
  euler_angles = {M_PI_2, 0., 0.};

  xp_yp_zp = euler_angle_transform::rotate(y_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), x_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, euler_angle_transform::rotate_back(xp_yp_zp, euler_angles).data(),
      y_axis.data(), epsilon);

  // Rotate y axis into z axis
  // Phi   = 0
  // Theta = -pi/2
  // Psi   = 0
  euler_angles = {0., -M_PI_2, 0.};

  xp_yp_zp = euler_angle_transform::rotate(y_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), z_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, euler_angle_transform::rotate_back(xp_yp_zp, euler_angles).data(),
      y_axis.data(), epsilon);

  // From here on, note the special role of the z axis in spherical coordinates
  // which leads to arbitrary values for phi.

  // Rotate z axis into x axis
  // Phi   = 0
  // Theta = pi/2
  // Psi   = pi/2
  euler_angles = {0., M_PI_2, M_PI_2};

  xp_yp_zp = euler_angle_transform::rotate(z_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), x_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, euler_angle_transform::rotate_back(xp_yp_zp, euler_angles).data(),
      z_axis.data(), epsilon);

  // Rotate z axis into y axis
  // Phi   = 0
  // Theta = pi/2
  // Psi   = 0
  euler_angles = {0., M_PI_2, 0.};

  xp_yp_zp = euler_angle_transform::rotate(z_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), y_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, euler_angle_transform::rotate_back(xp_yp_zp, euler_angles).data(),
      z_axis.data(), epsilon);

  // Rotate z axis into z axis (trivial)
  // Phi   = 0
  // Theta = 0
  // Psi   = 0
  euler_angles = {0., 0., 0.};

  xp_yp_zp = euler_angle_transform::rotate(z_axis, euler_angles);
  test_numerical_equality<double>(3, xp_yp_zp.data(), z_axis.data(), epsilon);
  test_numerical_equality<double>(
      3, euler_angle_transform::rotate_back(xp_yp_zp, euler_angles).data(),
      z_axis.data(), epsilon);

  // Test the get_theta_phi method which calculates the corresponding spherical
  // coordinates theta and phi for a given normalized Cartesian vector.
  array<double, 2> theta_phi;

  // Test the x, y, z, -x, -y, and -z axis.
  theta_phi = euler_angle_transform::get_theta_phi({1., 0., 0.});
  test_numerical_equality<double>(theta_phi[0], M_PI_2, epsilon);
  test_numerical_equality<double>(theta_phi[1], 0., epsilon);

  theta_phi = euler_angle_transform::get_theta_phi({-1., 0., 0.});
  test_numerical_equality<double>(theta_phi[0], M_PI_2, epsilon);
  test_numerical_equality<double>(theta_phi[1], M_PI, epsilon);

  theta_phi = euler_angle_transform::get_theta_phi({0., 1., 0.});
  test_numerical_equality<double>(theta_phi[0], M_PI_2, epsilon);
  test_numerical_equality<double>(theta_phi[1], M_PI_2, epsilon);

  theta_phi = euler_angle_transform::get_theta_phi({0., -1., 0.});
  test_numerical_equality<double>(theta_phi[0], M_PI_2, epsilon);
  test_numerical_equality<double>(theta_phi[1], 3. * M_PI_2, epsilon);

  theta_phi = euler_angle_transform::get_theta_phi({0., 0., 1.});
  test_numerical_equality<double>(theta_phi[0], 0., epsilon);
  test_numerical_equality<double>(theta_phi[1], 0., epsilon);

  theta_phi = euler_angle_transform::get_theta_phi({0., 0., -1.});
  test_numerical_equality<double>(theta_phi[0], M_PI, epsilon);
  test_numerical_equality<double>(theta_phi[1], 0., epsilon);

  // Test more vectors in the xy plane to see whether phi is set in the correct
  // quadrant.
  theta_phi = euler_angle_transform::get_theta_phi({1., 1., 0.});
  test_numerical_equality<double>(theta_phi[0], M_PI_2, epsilon);
  test_numerical_equality<double>(theta_phi[1], M_PI_4, epsilon);

  theta_phi = euler_angle_transform::get_theta_phi({-1., 1., 0.});
  test_numerical_equality<double>(theta_phi[0], M_PI_2, epsilon);
  test_numerical_equality<double>(theta_phi[1], 3. * M_PI_4, epsilon);

  theta_phi = euler_angle_transform::get_theta_phi({-1., -1., 0.});
  test_numerical_equality<double>(theta_phi[0], M_PI_2, epsilon);
  test_numerical_equality<double>(theta_phi[1], 5. * M_PI_4, epsilon);

  theta_phi = euler_angle_transform::get_theta_phi({1., -1., 0.});
  test_numerical_equality<double>(theta_phi[0], M_PI_2, epsilon);
  test_numerical_equality<double>(theta_phi[1], 7. * M_PI_4, epsilon);

  // Test back-and-forth conversion between Euler angles and coordinate axes.
  array<array<double, 3>, 3> transformed_axes;

  transformed_axes = euler_angle_transform::axes({M_PI_2, 0, 0});
  test_numerical_equality<double, 3>(
      transformed_axes,
      array<array<double, 3>, 3>{array<double, 3>{0, -1, 0},
                                 array<double, 3>{1, 0, 0},
                                 array<double, 3>{0, 0, 1}},
      epsilon);

  transformed_axes = euler_angle_transform::axes({0, M_PI_2, 0});
  test_numerical_equality<double, 3>(
      transformed_axes,
      array<array<double, 3>, 3>{array<double, 3>{1, 0, 0},
                                 array<double, 3>{0, 0, -1},
                                 array<double, 3>{0, 1, 0}},
      epsilon);

  transformed_axes = euler_angle_transform::axes({0, 0, M_PI_2});
  test_numerical_equality<double, 3>(
      transformed_axes,
      array<array<double, 3>, 3>{array<double, 3>{0, -1, 0},
                                 array<double, 3>{1, 0, 0},
                                 array<double, 3>{0, 0, 1}},
      epsilon);
}