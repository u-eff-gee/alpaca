/*
    This file is part of alpaca

    alpaca is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    alpaca is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with alpaca  If not, see <https://www.gnu.org/licenses/>.

    Copyright (C) 2021-2023 Udo Friman-Gayer
*/

#include <array>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include "EulerAngleRotation.hh"
#include "TestUtilities.hh"

void test_numerical_equality(const size_t n, gsl_vector *a, gsl_vector *b,
                             const double epsilon) {
  for (size_t i = 0; i < n; ++i) {
    test_numerical_equality<double>(gsl_vector_get(a, i), gsl_vector_get(b, i),
                                    epsilon);
  }
}

/**
 * \brief Test rotations of 3D vectors using Euler angles.
 *
 * Test by rotating the three canonical Cartesian axes into each other.
 */
int main() {

  const double epsilon = 1e-8;

  gsl_vector *x_axis = gsl_vector_alloc(3);
  gsl_vector_set_basis(x_axis, 0);
  gsl_vector *y_axis = gsl_vector_alloc(3);
  gsl_vector_set_basis(y_axis, 1);
  gsl_vector *z_axis = gsl_vector_alloc(3);
  gsl_vector_set_basis(z_axis, 2);

  // Rotate x axis into y axis
  gsl_vector *Phi_Theta_Psi_x_to_y = gsl_vector_alloc(3);
  gsl_vector_set_zero(Phi_Theta_Psi_x_to_y);
  gsl_vector_set(Phi_Theta_Psi_x_to_y, 0, -M_PI_2);

  gsl_vector *x_rotated_to_y = gsl_vector_alloc(3);
  euler_angle_transform::rotate(x_rotated_to_y, Phi_Theta_Psi_x_to_y, x_axis);
  test_numerical_equality(3, x_rotated_to_y, y_axis, epsilon);

  gsl_vector *x_rotated_to_y_reversed = gsl_vector_alloc(3);
  euler_angle_transform::rotate_back(x_rotated_to_y_reversed,
                                     Phi_Theta_Psi_x_to_y, x_rotated_to_y);
  test_numerical_equality(3, x_rotated_to_y_reversed, x_axis, epsilon);

  // Rotate x axis into z axis
  gsl_vector *Phi_Theta_Psi_x_to_z = gsl_vector_alloc(3);
  gsl_vector_set_zero(Phi_Theta_Psi_x_to_z);
  gsl_vector_set(Phi_Theta_Psi_x_to_z, 0, M_PI_2);
  gsl_vector_set(Phi_Theta_Psi_x_to_z, 1, M_PI_2);

  gsl_vector *x_rotated_to_z = gsl_vector_alloc(3);
  euler_angle_transform::rotate(x_rotated_to_z, Phi_Theta_Psi_x_to_z, x_axis);
  test_numerical_equality(3, x_rotated_to_z, z_axis, epsilon);

  gsl_vector *x_rotated_to_z_reversed = gsl_vector_alloc(3);
  euler_angle_transform::rotate_back(x_rotated_to_z_reversed,
                                     Phi_Theta_Psi_x_to_z, x_rotated_to_z);
  test_numerical_equality(3, x_rotated_to_z_reversed, x_axis, epsilon);

  // Rotate y axis into x axis
  gsl_vector *Phi_Theta_Psi_y_to_x = gsl_vector_alloc(3);
  gsl_vector_set_zero(Phi_Theta_Psi_y_to_x);
  gsl_vector_set(Phi_Theta_Psi_y_to_x, 0, M_PI_2);

  gsl_vector *y_rotated_to_x = gsl_vector_alloc(3);
  euler_angle_transform::rotate(y_rotated_to_x, Phi_Theta_Psi_y_to_x, y_axis);
  test_numerical_equality(3, y_rotated_to_x, x_axis, epsilon);

  gsl_vector *y_rotated_to_x_reversed = gsl_vector_alloc(3);
  euler_angle_transform::rotate_back(y_rotated_to_x_reversed,
                                     Phi_Theta_Psi_y_to_x, y_rotated_to_x);
  test_numerical_equality(3, y_rotated_to_x_reversed, y_axis, epsilon);

  // Rotate y axis into z axis
  gsl_vector *Phi_Theta_Psi_y_to_z = gsl_vector_alloc(3);
  gsl_vector_set_zero(Phi_Theta_Psi_y_to_z);
  gsl_vector_set(Phi_Theta_Psi_y_to_z, 1, -M_PI_2);

  gsl_vector *y_rotated_to_z = gsl_vector_alloc(3);
  euler_angle_transform::rotate(y_rotated_to_z, Phi_Theta_Psi_y_to_z, y_axis);
  test_numerical_equality(3, y_rotated_to_z, z_axis, epsilon);

  gsl_vector *y_rotated_to_z_reversed = gsl_vector_alloc(3);
  euler_angle_transform::rotate_back(y_rotated_to_z_reversed,
                                     Phi_Theta_Psi_y_to_z, y_rotated_to_z);
  test_numerical_equality(3, y_rotated_to_z_reversed, y_axis, epsilon);

  // Rotate z axis into x axis
  gsl_vector *Phi_Theta_Psi_z_to_x = gsl_vector_alloc(3);
  gsl_vector_set_zero(Phi_Theta_Psi_z_to_x);
  gsl_vector_set(Phi_Theta_Psi_z_to_x, 1, M_PI_2);
  gsl_vector_set(Phi_Theta_Psi_z_to_x, 2, M_PI_2);

  gsl_vector *z_rotated_to_x = gsl_vector_alloc(3);
  euler_angle_transform::rotate(z_rotated_to_x, Phi_Theta_Psi_z_to_x, z_axis);
  test_numerical_equality(3, z_rotated_to_x, x_axis, epsilon);

  gsl_vector *z_rotated_to_x_reversed = gsl_vector_alloc(3);
  euler_angle_transform::rotate_back(z_rotated_to_x_reversed,
                                     Phi_Theta_Psi_z_to_x, z_rotated_to_x);
  test_numerical_equality(3, z_rotated_to_x_reversed, z_axis, epsilon);

  // Rotate z axis into y axis
  gsl_vector *Phi_Theta_Psi_z_to_y = gsl_vector_alloc(3);
  gsl_vector_set_zero(Phi_Theta_Psi_z_to_y);
  gsl_vector_set(Phi_Theta_Psi_z_to_y, 1, M_PI_2);

  gsl_vector *z_rotated_to_y = gsl_vector_alloc(3);
  euler_angle_transform::rotate(z_rotated_to_y, Phi_Theta_Psi_z_to_y, z_axis);
  test_numerical_equality(3, z_rotated_to_y, y_axis, epsilon);

  gsl_vector *z_rotated_to_y_reversed = gsl_vector_alloc(3);
  euler_angle_transform::rotate_back(z_rotated_to_y_reversed,
                                     Phi_Theta_Psi_z_to_y, z_rotated_to_y);
  test_numerical_equality(3, z_rotated_to_y_reversed, z_axis, epsilon);
}