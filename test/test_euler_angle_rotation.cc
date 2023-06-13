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

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/TestUtilities.hh"

using std::exception;
using std::runtime_error;
using std::string;
using std::stringstream;
using std::vector;

using alpaca::EulerAngles;
using alpaca::test_numerical_equality;
namespace euler_angle_transform = alpaca::euler_angle_transform;

void test_numerical_equality(gsl_vector *a, gsl_vector *b,
                             const double epsilon) {
  for (size_t i = 0; i < 3; ++i) {
    try {
      test_numerical_equality<double>(gsl_vector_get(a, i),
                                      gsl_vector_get(b, i), epsilon);
    } catch (const exception &e) {
      stringstream error_message;
      error_message << e.what();
      error_message << "For index " << i << " of vectors"
                    << "\n";
      error_message << "(" << std::scientific << gsl_vector_get(a, 0) << ", "
                    << gsl_vector_get(a, 1) << ", " << gsl_vector_get(a, 2)
                    << ")\n";
      error_message << "(" << std::scientific << gsl_vector_get(b, 0) << ", "
                    << gsl_vector_get(b, 1) << ", " << gsl_vector_get(b, 2)
                    << ")\n";
      throw runtime_error(error_message.str());
    }
  }
}

struct RotationTest {
  string name;
  EulerAngles initial_axis;
  EulerAngles target_axis;
  EulerAngles Phi_Theta_Psi;
  EulerAngles Phi_Theta_Psi_alternative;
};

const vector<RotationTest> rotation_tests{
    // If Phi_Theta_Psi_alternative is initialized with zeros only, that means
    // it won't be used.
    RotationTest{"x_to_y", {1, 0, 0}, {0, 1, 0}, {M_PI_2, 0, 0}, {0., 0., 0.}},
    RotationTest{
        "x_to_z", {1, 0, 0}, {0, 0, 1}, {M_PI_2, M_PI_2, 0}, {0., 0., 0.}},
    RotationTest{"y_to_x", {0, 1, 0}, {1, 0, 0}, {-M_PI_2, 0, 0}, {0., 0., 0.}},
    RotationTest{
        "y_to_z", {0, 1, 0}, {0, 0, 1}, {M_PI, -M_PI_2, M_PI}, {0, M_PI_2, 0}},
    RotationTest{"z_to_x",
                 {0, 0, 1},
                 {1, 0, 0},
                 {0, -M_PI_2, -M_PI_2},
                 {-M_PI, M_PI_2, M_PI_2}},
    RotationTest{"z_to_y",
                 {0, 0, 1},
                 {0, 1, 0},
                 {0, -M_PI_2, 0},
                 {-M_PI, M_PI_2, -M_PI}},
    RotationTest{"x_to_x", {1, 0, 0}, {1, 0, 0}, {0, 0, 0}, {0, 0, 0}},
};

/**
 * \brief Test rotations of 3D vectors using Euler angles.
 *
 * This program employs user-defined sets of Euler angles (1) to transform the
 * primitive x, y, and z axes into each other (and one case where no
 * transformation is necessary). After that the backward transformation is used
 * to restore the initial vector (2). In a third step (3), the rotation matrix
 * for the transformation in (1) is explicitly calculated, and it is
 * demonstrated how to infer the Euler angles from it. Since the reverse
 * engineering of Euler angles is not unique, the user-defined Euler angles may
 * provide an alternative set of angles (knowing what will be the result of
 * euler_angle_transform::angles) to achieve the same transformation. This
 * ambiguity is handled by a try-catch clause (4). At least, it is verified that
 * the original and alternative sets of Euler angles indeed describe the same
 * transformation by applying it too all three primitive axes (5).
 */
int main() {

  const double epsilon = 1e-8;

  gsl_vector *x_axis = gsl_vector_alloc(3);
  gsl_vector_set_basis(x_axis, 0);
  gsl_vector *x_axis_transformed = gsl_vector_alloc(3);
  gsl_vector *x_axis_transformed_from_reconstructed_matrix =
      gsl_vector_alloc(3);
  gsl_vector *y_axis = gsl_vector_alloc(3);
  gsl_vector_set_basis(y_axis, 1);
  gsl_vector *y_axis_transformed = gsl_vector_alloc(3);
  gsl_vector *y_axis_transformed_from_reconstructed_matrix =
      gsl_vector_alloc(3);
  gsl_vector *z_axis = gsl_vector_alloc(3);
  gsl_vector_set_basis(z_axis, 2);
  gsl_vector *z_axis_transformed = gsl_vector_alloc(3);
  gsl_vector *z_axis_transformed_from_reconstructed_matrix =
      gsl_vector_alloc(3);

  gsl_vector *initial_axis = gsl_vector_alloc(3);
  gsl_vector *Phi_Theta_Psi = gsl_vector_alloc(3);
  gsl_vector *initial_axis_transformed = gsl_vector_alloc(3);
  gsl_vector *initial_axis_reconstructed = gsl_vector_alloc(3);
  gsl_vector *target_axis = gsl_vector_alloc(3);
  gsl_matrix *transformation_matrix = gsl_matrix_alloc(3, 3);
  gsl_vector *Phi_Theta_Psi_reconstructed = gsl_vector_alloc(3);
  gsl_vector *Phi_Theta_Psi_alternative = gsl_vector_alloc(3);

  for (auto test : rotation_tests) {
    gsl_vector_set(initial_axis, 0, test.initial_axis[0]);
    gsl_vector_set(initial_axis, 1, test.initial_axis[1]);
    gsl_vector_set(initial_axis, 2, test.initial_axis[2]);

    gsl_vector_set(target_axis, 0, test.target_axis[0]);
    gsl_vector_set(target_axis, 1, test.target_axis[1]);
    gsl_vector_set(target_axis, 2, test.target_axis[2]);

    gsl_vector_set(Phi_Theta_Psi, 0, test.Phi_Theta_Psi[0]);
    gsl_vector_set(Phi_Theta_Psi, 1, test.Phi_Theta_Psi[1]);
    gsl_vector_set(Phi_Theta_Psi, 2, test.Phi_Theta_Psi[2]);

    // (1) Verify that the given transformation turns the given initial axis
    // into the given target axis.
    euler_angle_transform::rotate(initial_axis_transformed, Phi_Theta_Psi,
                                  initial_axis);
    test_numerical_equality(initial_axis_transformed, target_axis, epsilon);

    // (2) Verify that the euler_angle_transform::rotate_back function does
    // exactly the opposite of euler_angle_transform::rotate.
    euler_angle_transform::rotate_back(initial_axis_reconstructed,
                                       Phi_Theta_Psi, initial_axis_transformed);
    test_numerical_equality(initial_axis_reconstructed, initial_axis, epsilon);

    // (3) Obtain forward-transformation matrix and reconstruct the Euler angles
    // from it.
    euler_angle_transform::rotation_matrix(transformation_matrix,
                                           Phi_Theta_Psi);
    euler_angle_transform::angles(Phi_Theta_Psi_reconstructed,
                                  transformation_matrix);
    try {
      // (4) Check whether the reconstructed Euler angles are the same as the
      // original ones. ...
      test_numerical_equality(Phi_Theta_Psi_reconstructed, Phi_Theta_Psi,
                              epsilon);
    } catch (exception &e) {
      // (4) ... Sometimes, this is not the case, because the transformations
      // are ambiguous. For this reason, the test cases contain an alternative
      // set of Euler angles. The alternative set was found -by trial and error-
      // to be the one that is really returned by euler_angle_transform::angles.
      gsl_vector_set(Phi_Theta_Psi_alternative, 0,
                     test.Phi_Theta_Psi_alternative[0]);
      gsl_vector_set(Phi_Theta_Psi_alternative, 1,
                     test.Phi_Theta_Psi_alternative[1]);
      gsl_vector_set(Phi_Theta_Psi_alternative, 2,
                     test.Phi_Theta_Psi_alternative[2]);
      test_numerical_equality(Phi_Theta_Psi_reconstructed,
                              Phi_Theta_Psi_alternative, epsilon);
    }

    // (5) Verify that the alternative set of Euler angles and the original one
    // indeed perform the same transformation by applying it to all three
    // primitive axes.
    euler_angle_transform::rotate(x_axis_transformed, Phi_Theta_Psi, x_axis);
    euler_angle_transform::rotate(x_axis_transformed_from_reconstructed_matrix,
                                  Phi_Theta_Psi_reconstructed, x_axis);
    test_numerical_equality(x_axis_transformed_from_reconstructed_matrix,
                            x_axis_transformed, epsilon);
    euler_angle_transform::rotate(y_axis_transformed, Phi_Theta_Psi, y_axis);
    euler_angle_transform::rotate(y_axis_transformed_from_reconstructed_matrix,
                                  Phi_Theta_Psi_reconstructed, y_axis);
    test_numerical_equality(y_axis_transformed_from_reconstructed_matrix,
                            y_axis_transformed, epsilon);
    euler_angle_transform::rotate(z_axis_transformed, Phi_Theta_Psi, z_axis);
    euler_angle_transform::rotate(z_axis_transformed_from_reconstructed_matrix,
                                  Phi_Theta_Psi_reconstructed, z_axis);
    test_numerical_equality(z_axis_transformed_from_reconstructed_matrix,
                            z_axis_transformed, epsilon);
  }

  gsl_vector_free(x_axis);
  gsl_vector_free(x_axis_transformed);
  gsl_vector_free(x_axis_transformed_from_reconstructed_matrix);
  gsl_vector_free(y_axis);
  gsl_vector_free(y_axis_transformed);
  gsl_vector_free(y_axis_transformed_from_reconstructed_matrix);
  gsl_vector_free(z_axis);
  gsl_vector_free(z_axis_transformed);
  gsl_vector_free(z_axis_transformed_from_reconstructed_matrix);

  gsl_vector_free(initial_axis);
  gsl_vector_free(Phi_Theta_Psi);
  gsl_vector_free(initial_axis_transformed);
  gsl_vector_free(initial_axis_reconstructed);
  gsl_vector_free(target_axis);
  gsl_matrix_free(transformation_matrix);
  gsl_vector_free(Phi_Theta_Psi_reconstructed);
  gsl_vector_free(Phi_Theta_Psi_alternative);
}
