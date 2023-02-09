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

#include "EulerAngleRotation.hh"

/**
 * @brief Abstract class for sampling a direction vector in spherical
 * coordinates.
 *
 * It is assumed that derived classes sample a direction vector with the
 * spherical coordinates \f$\theta\f$ (polar angle) and \f$\varphi\f$ (azimuthal
 * angle) using an algortihm that may require a number of tries before a valid
 * vector is found [e.g. rejection sampling (see, e.g. Sec. 2.3
 * in Ref. \cite RobertCasella1999 and the *RejectionSampler classes
 * for more information].
 *
 * This class defines the abstract interface for the samplers, including a
 * method to report the number of tries until a valid vector was found and
 * another one to rotate the sampled vector into a different coordinate system.
 */
class DirectionSampler {
public:
  DirectionSampler() : euler_angle_rotation(EulerAngleRotation()) {}

  /**
   * \brief Sample a random vector and record the
   * number of tries.
   *
   * \return std::pair which contains \f$N\f$, the number of tries that were
   * needed to find a valid vector, and the accepted vector \f$\left(
   * \theta_\mathrm{rand}, \varphi_\mathrm{rand}\right)\f$. Returns a std::pair
   * of \f$N_\mathrm{max}\f$ and \f$\left(0, 0 \right)\f$, if the maximum number
   * of trials \f$N_\mathrm{max}\f$ is reached by the algorithm and no random
   * vector was accepted.
   */
  virtual pair<unsigned int, array<double, 2>> sample() = 0;
  /**
   * \brief Sample a random vector.
   *
   * \return Accepted vector \f$\left( \theta_\mathrm{rand},
   * \varphi_\mathrm{rand}\right)\f$. Returns \f$\left( 0, 0 \right)\f$ if the
   * maximum number of trials \f$N_\mathrm{max}\f$ is reached by the algorithm
   * and no random vector was accepted.
   */
  array<double, 2> operator()();

  /**
   * \brief Sample a random vector from an arbitrarily rotated probability
   * distribution.
   *
   * This function allows to rotate the probability distribution using Euler
   * angles in the x convention. In the present code, this is an important
   * feature in gamma-ray cascades, where the direction of emission/propagation
   * and the polarization axis of the initial photon define the reference frame
   * for the emission of a subsequent photon. [Usually the analytical
   * expressions for these angular correlations are implemented in a fixed
   * reference frame (here, the direction of emission/propagation is along the
   * positive \f$z\f$ axis, and the polarization along the \f$x\f$ axis) for
   * simplicity.]
   *
   * \param euler_angles Euler angles \f$\Phi\f$, \f$\Theta\f$, and \f$\Psi\f$
   * in radians which define an arbitrary rotation in 3D space in the x
   * convention.
   *
   * \return Accepted vector \f$\left( \theta_\mathrm{rand},
   * \varphi_\mathrm{rand}\right)\f$. Returns \f$\left( 0, 0 \right)\f$ if the
   * maximum number of trials \f$N_\mathrm{max}\f$ is reached by the algorithm
   * and no random vector was accepted.
   */
  array<double, 2> operator()(const array<double, 3> euler_angles);

protected:
  const EulerAngleRotation
      euler_angle_rotation; /**< Instance of the EulerAngleRotation class */
};