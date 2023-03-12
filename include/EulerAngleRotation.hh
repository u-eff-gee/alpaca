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
 * For the representation of 2D and 3D vectors, this class uses the std::array
 * container class.
 *
 */
namespace euler_angle_transform {

/**
 * \brief Check if all three Euler angles are zero.
 *
 * This function is used in the code to decide whether a calculation needs to
 * be performed or whether the input value can simply be returned.
 *
 * \param Phi_Theta_Psi Euler angles in radians
 *
 * \return true, if all Euler angles are zero, else false.
 */
inline bool no_rotation_required(const array<double, 3> Phi_Theta_Psi) {
  if (Phi_Theta_Psi[0] == 0. && Phi_Theta_Psi[1] == 0. &&
      Phi_Theta_Psi[2] == 0.) {
    return true;
  }
  return false;
};

/**
 * \brief Calculate rotation matrix for the three Euler angles.
 *
 * \param Phi_Theta_Psi Euler angles in radians
 *
 * \return \f$3 \times 3\f$ matrix \f$A\f$
 */
inline array<array<double, 3>, 3>
rotation_matrix(const array<double, 3> Phi_Theta_Psi) {

  const double cos_phi{cos(Phi_Theta_Psi[0])}, sin_phi{sin(Phi_Theta_Psi[0])},
      cos_the{cos(Phi_Theta_Psi[1])}, sin_the{sin(Phi_Theta_Psi[1])},
      cos_psi{cos(Phi_Theta_Psi[2])}, sin_psi{sin(Phi_Theta_Psi[2])};

  return array<array<double, 3>, 3>{
      array<double, 3>{cos_psi * cos_phi - cos_the * sin_phi * sin_psi,
                       cos_psi * sin_phi + cos_the * cos_phi * sin_psi,
                       sin_psi * sin_the},
      array<double, 3>{-sin_psi * cos_phi - cos_the * sin_phi * cos_psi,
                       -sin_psi * sin_phi + cos_the * cos_phi * cos_psi,
                       cos_psi * sin_the},
      array<double, 3>{sin_the * sin_phi, -sin_the * cos_phi, cos_the}};
};

/**
 * \brief Convert spherical to Cartesian coordinates.
 *
 * Given spherical coordinates \f$\theta\f$ and \f$\varphi\f$, this function
 * calculates the corresponding normalized Cartesian three-component vector.
 *
 * \param theta_phi Spherical coordinates \f$\theta\f$ and \f$\varphi\f$ in
 * radians.
 *
 * \return Normalized Cartesian vector.
 */
inline array<double, 3> get_x_y_z_norm(const array<double, 2> theta_phi) {

  double cos_the{cos(theta_phi[0])}, sin_the{sin(theta_phi[0])},
      cos_phi{cos(theta_phi[1])}, sin_phi{sin(theta_phi[1])};

  return array<double, 3>{sin_the * cos_phi, sin_the * sin_phi, cos_the};
};

/**
 * \brief Convert Cartesian to spherical coordinates.
 *
 * Given a normalized Cartesian vector with three coordinates \f$x\f$,
 * \f$y\f$, and \f$z\f$, this function calculates the corresponding angles
 * \f$\theta\f$ and \f$\varphi\f$ in spherical coordinates.
 *
 * At the moment, the input vector is not tested for normalization.
 *
 * \param x_y_z_norm Normalized Cartesian vector.
 *
 * \return Spherical coordinates \f$\theta\f$ and \f$\varphi\f$ in radians.
 */
inline array<double, 2> get_theta_phi(const array<double, 3> x_y_z_norm) {

  return {acos(x_y_z_norm[2]),
          fmod(atan2(x_y_z_norm[1], x_y_z_norm[0]) + 2. * M_PI, 2. * M_PI)};
};

/**
 * \brief Rotate a 3D vector.
 *
 * Vector to rotate is given in Cartesian coordinates (and does not have to be
 * normalized), rotation is represented as a rotation matrix.
 *
 * \param x_y_z \f$v\f$, 3D vector
 * \param A \f$A\$, 3 x 3 rotation matrix.
 *
 * \return \f$v^\prime\f$, 3D vector
 */
inline array<double, 3> rotate(const array<double, 3> x_y_z,
                        const array<array<double, 3>, 3> A) {
  return array<double, 3>{
      A[0][0] * x_y_z[0] + A[0][1] * x_y_z[1] + A[0][2] * x_y_z[2],
      A[1][0] * x_y_z[0] + A[1][1] * x_y_z[1] + A[1][2] * x_y_z[2],
      A[2][0] * x_y_z[0] + A[2][1] * x_y_z[1] + A[2][2] * x_y_z[2]};
};

/**
 * \brief Rotate a 3D vector.
 *
 * Vector to rotate is given in spherical coordinates, rotation is represented
 * as set of Euler angles.
 *
 * \param theta_phi \f$v\f$, spherical coordinates \f$\theta\f$ and
 * \f$\varphi\f$ in radians.
 * \param A \f$A\$, 3 x 3 rotation matrix.
 *
 * \return \f$v^\prime\f$, spherical coordinates \f$\theta\f$ and
 * \f$\varphi\f$ in radians.
 */
inline array<double, 2> rotate(const array<double, 2> theta_phi,
                        const array<array<double, 3>, 3> A) {
  return get_theta_phi(rotate(get_x_y_z_norm(theta_phi), A));
};

/**
 * \brief Rotate a 3D vector.
 *
 * Vector to rotate is given in Cartesian coordinates (and does not have to be
 * normalized), rotation is represented as set of Euler angles.
 *
 * \param x_y_z \f$v\f$, 3D vector
 * \param phi_theta_psi Euler angles in radians
 *
 * \return \f$v^\prime\f$, 3D vector
 */
inline array<double, 3> rotate(const array<double, 3> x_y_z,
                        const array<double, 3> Phi_Theta_Psi) {

  if (euler_angle_transform::no_rotation_required(Phi_Theta_Psi)) {
    return x_y_z;
  }

  const array<array<double, 3>, 3> A = rotation_matrix(Phi_Theta_Psi);

  return rotate(x_y_z, A);
};

/**
 * \brief Rotate a 3D vector.
 *
 * Vector to rotate is given in spherical coordinates, rotation is represented
 * as a rotation matrix..
 *
 * \param theta_phi \f$v\f$, spherical coordinates \f$\theta\f$ and
 * \f$\varphi\f$ in radians. \param Phi_Theta_Psi Euler angles in radians
 *
 * \return \f$v^\prime\f$, spherical coordinates \f$\theta\f$ and
 * \f$\varphi\f$ in radians.
 */
inline array<double, 2> rotate(const array<double, 2> theta_phi,
                        const array<double, 3> Phi_Theta_Psi) {
  if (euler_angle_transform::no_rotation_required(Phi_Theta_Psi)) {
    return theta_phi;
  }

  return get_theta_phi(rotate(get_x_y_z_norm(theta_phi), Phi_Theta_Psi));
};

/**
 * \brief Rotate a 3D vector back.
 *
 * Performs the same action on a 3D vector which would be performed by
 * EulerAngleRotation::rotate() if the angles \f$\Phi\f$ and \f$\Psi\f$ are
 * switched, and the negative value of each angle is used.
 *
 * \param x_y_z \f$v^\prime\f$, 3D vector
 * \param Phi_Theta_Psi Euler angles in radians
 *
 * \return \f$v\f$, 3D vector
 */
inline array<double, 3> rotate_back(const array<double, 3> xp_yp_zp,
                             const array<double, 3> Phi_Theta_Psi) {
  if (euler_angle_transform::no_rotation_required(Phi_Theta_Psi)) {
    return xp_yp_zp;
  }

  return rotate(xp_yp_zp, array<double, 3>{-Phi_Theta_Psi[2], -Phi_Theta_Psi[1],
                                           -Phi_Theta_Psi[0]});
};

/**
 * \brief Rotate a 3D vector back.
 *
 * See also the implementation of rotate_back for Cartesian vectors.
 *
 * \param theta_phi \f$v^\prime\f$, spherical coordinates \f$\theta\f$ and
 * \f$\varphi\f$ in radians. \param Phi_Theta_Psi Euler angles in radians
 *
 * \return \f$v\f$, spherical coordinates \f$\theta\f$ and \f$\varphi\f$ in
 * radians.
 */
inline array<double, 2> rotate_back(const array<double, 2> thetap_phip,
                             const array<double, 3> Phi_Theta_Psi) {
  if (euler_angle_transform::no_rotation_required(Phi_Theta_Psi)) {
    return thetap_phip;
  }

  return get_theta_phi(rotate_back(get_x_y_z_norm(thetap_phip), Phi_Theta_Psi));
};
}; // namespace euler_angle_transform
