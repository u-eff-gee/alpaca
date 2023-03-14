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

#include <array>

using std::array;

#include <cmath>

#include <gsl/gsl_blas.h>

/**
 * \brief Functions to perform arbitrary rotations of 3D vectors using Euler
 * angles.
 *
 * Any rotation in three-dimensional space can be expressed in terms of three
 * Euler angles, which are denoted as \f$\Phi\f$, \f$\Theta\f$, and \f$\Psi\f$
 * here \cite Weisstein2020. Here, an arbitrary orientation of a coordinate
 * system with the axes \f$x\f$, \f$y\f$, and \f$z\f$ into a coordinate system
 * with the new axes \f$x^\prime\f$, \f$y^\prime\f$, and \f$z^\prime\f$ is
 * achieved by the 'zxz' scheme or the 'x-convention': First, the original
 * vector \f$v\f$ in the \f$xyz\f$ system is rotated around the original \f$z\f$
 * axis by an angle \f$\Phi\f$. The rotation matrix is given by:
 *
 * \f[
 *      D\left( \Phi \right) = \begin{pmatrix}
 *           \cos \left( \Phi \right) &  \sin \left( \Phi \right) & 0 \\
 *          -\sin \left( \Phi \right) &  \cos \left( \Phi \right) & 0 \\
 *           0                        &  0                        & 1 \\
 *      \end{pmatrix}
 * \f]
 *
 * Then, the resulting vector is rotated around the original \f$x\f$ axis by an
 * angle \f$\Theta\f$:
 *
 * \f[
 *      C\left( \Theta \right) = \begin{pmatrix}
 *           1                        &  0                          & 0 \\
 *           0                        &  \cos \left( \Theta \right) &  \sin
 * \left( \Theta \right) \\
 *           0                        & -\sin \left( \Theta \right) &  \cos
 * \left( \Theta \right) \\ \end{pmatrix} \f]
 *
 * At last, the vector is rotated around the original \f$z\f$ axis by an angle
 * \f$\Psi\f$.
 *
 * \f[
 *      B\left( \Psi \right) = \begin{pmatrix}
 *           \cos \left( \Psi \right) &  \sin \left( \Psi \right) & 0 \\
 *          -\sin \left( \Psi \right) &  \cos \left( \Psi \right) & 0 \\
 *           0                        &  0                        & 1 \\
 *      \end{pmatrix}
 * \f]
 *
 * The matrices are named in the same way as Eq. (1) in Ref. \cite
 * Weisstein2020, but they have explicit angle arguments. The action of the
 * three rotation matrices on the vector \f$v\f$ results in a new vector
 * \f$v^\prime\f$:
 *
 * \f[
 *      v^\prime = \underbrace{B \left( \Psi \right) C \left( \Theta \right) D
 * \left( \Phi \right)}_{\equiv A \left( \Phi, \Theta, \Psi \right)} v. \f]
 *
 * In the equation above, the symbol \f$A\f$ for the total rotation matrix has
 * been introduced. In order to obtain \f$v\f$ from \f$v^\prime\f$, the
 * rotations can be reversed:
 *
 * \f[
 *      v = \underbrace{D \left( -\Phi \right) C \left( -\Theta \right) B \left(
 * -\Psi \right)}_{\equiv A^{-1} \left( \Phi, \Theta, \Psi \right)} v^\prime.
 * \f]
 *
 * This is equivalent to calculating the inverse matrix \f$A^{-1}\f$ of \f$A\f$.
 *
 * For the present application, it is instructive to consider the action of the
 * three matrices on a unit vector along the positive \f$z\f$ axis, since the
 * default orientation of all angular correlations is along the \f$z\f$ axis (in
 * other words: the first photon is supposed to be emitted/propagating along the
 * \f$z\f$ axis):
 *
 * \f[
 *      B \left( \Psi \right) C \left( \Theta \right) D \left( \Phi \right)
 * \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} = B \left( \Psi \right) C \left(
 * \Theta \right) \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix}
 * \sin \left( \Theta \right) \sin \left( \Psi \right) \\ \sin \left( \Theta
 * \right) \cos \left( \Psi \right) \\ \cos \left( \Theta \right) \end{pmatrix}.
 * \f]
 *
 * For the convention used in this code, an arbitrary vector in spherical
 * coordinates is identified by its polar angle \f$\theta\f$ and its azimuthal
 * angle \f$\varphi\f$:
 *
 * \f[
 *      \begin{pmatrix} \sin \left( \theta \right) \cos \left( \varphi \right)
 * \\ \sin \left( \theta \right) \sin \left( \varphi \right) \\ \cos \left(
 * \theta \right) \end{pmatrix}. \f]
 *
 * By comparing the two equations above, it can be seen that a given vector in
 * spherical coordinates can be obtained by applying the matrix \f$B C\f$ to a
 * unit vector along the positive \f$z\f$ axis with the following choice of the
 * Euler angles:
 *
 * \f[
 *      \Theta = \theta
 * \f]
 * \f[
 *      \Psi = -\varphi + \frac{\pi}{2}.
 * \f]
 *
 * For the representation and manipulation of vectors and matrices, this class
 * uses GSL \cite Galassi2009 and its BLAS \cite Lawson1979 \cite Dongarra1988
 * \cite Dongarra1990 interface.
 *
 */
namespace euler_angle_transform {

inline void rotation_matrix(gsl_matrix *A, gsl_vector *Phi_Theta_Psi) {

  const double cos_phi{cos(Phi_Theta_Psi->data[0])},
      sin_phi{sin(Phi_Theta_Psi->data[0])},
      cos_the{cos(Phi_Theta_Psi->data[1])},
      sin_the{sin(Phi_Theta_Psi->data[1])},
      cos_psi{cos(Phi_Theta_Psi->data[2])},
      sin_psi{sin(Phi_Theta_Psi->data[2])};

  gsl_matrix_set(A, 0, 0, cos_psi * cos_phi - cos_the * sin_phi * sin_psi);
  gsl_matrix_set(A, 0, 1, cos_psi * sin_phi + cos_the * cos_phi * sin_psi);
  gsl_matrix_set(A, 0, 2, sin_psi * sin_the);
  gsl_matrix_set(A, 1, 0, -sin_psi * cos_phi - cos_the * sin_phi * cos_psi);
  gsl_matrix_set(A, 1, 1, -sin_psi * sin_phi + cos_the * cos_phi * cos_psi);
  gsl_matrix_set(A, 1, 2, cos_psi * sin_the);
  gsl_matrix_set(A, 2, 0, sin_the * sin_phi);
  gsl_matrix_set(A, 2, 1, -sin_the * cos_phi);
  gsl_matrix_set(A, 2, 2, cos_the);
};

inline gsl_vector angles(const gsl_matrix A) {
  double v[] = {
      atan2(gsl_matrix_get(&A, 2, 0), gsl_matrix_get(&A, 2, 1)),
      acos(gsl_matrix_get(&A, 2, 2)),
      atan2(gsl_matrix_get(&A, 0, 2), -gsl_matrix_get(&A, 1, 2)),
  };

  return gsl_vector_view_array(v, 3).vector;
}

inline void rotate(gsl_vector *xp_yp_zp, gsl_vector *Phi_Theta_Psi,
                   const gsl_vector *x_y_z) {

  gsl_matrix *A = gsl_matrix_alloc(3, 3);
  rotation_matrix(A, Phi_Theta_Psi);

  gsl_blas_dgemv(CblasNoTrans, 1., A, x_y_z, 0., xp_yp_zp);
};

inline void rotate_back(gsl_vector *x_y_z, gsl_vector *Phi_Theta_Psi,
                        const gsl_vector *xp_yp_zp) {

  gsl_vector_reverse(Phi_Theta_Psi);
  gsl_vector_scale(Phi_Theta_Psi, -1);

  rotate(x_y_z, Phi_Theta_Psi, xp_yp_zp);
};

/**
 * @brief Euler angles for rotating the z axis into an arbitrary vector in
 * spherical coordinates.
 *
 * Given a polar angle \f$\theta\f$ and and azimuthal angle \f$\varphi\f$ in
 * spherical coordinates, this function returns three Euler angles for a
 * rotation which turns the z axis (\f$\theta = 0\f$) into the given vector.
 * Since there are three Euler angles, but only two are needed for the
 * transformation, the third angle can be set by the user. That third angle is
 * supposed to be the one associated with the first rotation about the z axis,
 * \f$\Phi\f$, since it has no effect on the z axis.
 *
 * @param theta_phi Polar- and azimuthal angle in spherical coordinates in
 * radians.
 * @param Phi Euler angle for the first rotation around the z axis in radians.
 * @return array<double, 3> One possible set of Euler angles to rotate the z
 * axis into the given vector in spherical coordinates.
 */
inline array<double, 3> from_spherical(const array<double, 2> theta_phi,
                                       const double Phi = 0.) {
  return {Phi, theta_phi[0], M_PI_2 - theta_phi[1]};
}

/**
 * @brief New orientation of z-axis (z'') in spherical coordinates after
 * Euler-angle rotation.
 *
 * @param Phi_Theta_Psi Euler angles in radians.
 * radians.
 * @return array<double, 2> Polar and azimuthal angle of the z axis after
 * rotation by three Euler angles in radians.
 */
inline array<double, 2> to_spherical(const array<double, 3> Phi_Theta_Psi) {
  return {Phi_Theta_Psi[1], M_PI_2 - Phi_Theta_Psi[2]};
}

}; // namespace euler_angle_transform
