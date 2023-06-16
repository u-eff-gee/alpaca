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

#include <functional>
#include <random>
#include <utility>

using std::function;
using std::mt19937;
using std::pair;
using std::uniform_real_distribution;

#include "alpaca/ReferenceFrameSampler.hh"

namespace alpaca {

/**
 * \brief Sample from a probability distribution in spherical coordinates using
 * rejection sampling.
 *
 * Given a probability distribution
 *
 * \f[
 * W\left( \Omega \right) \mathrm{d}\Omega = W \left( \theta, \varphi\right)
 * \sin \left( \theta \right) \mathrm{d} \theta \mathrm{d} \varphi, \f]
 *
 * which determines the probability of finding a vector in the solid angle
 * element \f$\mathrm{d}\Omega\f$ (expressed as \f$\sin \left( \theta \right)
 * \mathrm{d} \theta \mathrm{d} \varphi\f$ in spherical coordinates, with the
 * polar angle \f$\theta\f$ and the aximuthal angle \f$\varphi\f$) this class
 * samples random Euler angles \f$\left( \Phi, \Theta, \Psi \right)\f$
 * distributed according to \f$W\f$.
 * The action of Euler-angle rotations by the aforementioned set of three angles
 * on the canonical z axis \f$\left( 0, 0, 1\right)\f$ will give vectors that
 * are distributed according to \f$W\f$.
 * When sampling spherical coordinates like this, it is clear that the value of
 * the first Euler angle \f$\Phi\f$ is irrelevant. It describes a rotation
 * around the original z axis, which has no effect on the z axis. Therefore,
 * SphereRejectionSampler initializes \f$\Phi\f$ with an random value in the
 * range \f$\left[ \right. 0, 2 \pi \left. \right)\f$.
 * In the language of polarization-direction correlations, this means that only
 * the direction of a particle is measured, but not the polarization, resulting
 * in a cylindrical symmetry around the z axis.
 *
 * The procedure above is equivalent to sampling a random rotation in 3D, which
 * is described, for example, in \cite Arvo1991 (which is the base for the more
 * easily accessible Ref. \cite Becker2012):
 *
 * "To generate uniformly distributed random rotations of a unit sphere, first
 * perform a random rotation about the vertical axis, then rotate the north pole
 * to a random position."
 *
 * The sampling algorithm used here is 'rejection sampling' (see, e.g. Sec. 2.3
 * in Ref. \cite RobertCasella1999). It requires an upper limit
 * \f$W_\mathrm{max}\f$ for the maximum value of the distribution:
 *
 * \f[
 *      W_\mathrm{max} \geq \mathrm{max}_{\theta \in \left[ 0, \pi \right],
 * \varphi \in \left[ 0, 2\pi \right]} W \left( \theta, \varphi \right). \f]
 *
 * It starts by sampling a point \f$\left( \theta_\mathrm{rand},
 * \varphi_\mathrm{rand} \right)\f$ from a uniform distribution on a sphere
 * surface (see, e.g. Ref. \cite Weisstein20202):
 *
 * \f[
 *      \theta_\mathrm{rand} = \mathrm{arccos} \left( 2u - 1 \right)
 * \f]
 * \f[
 *      \varphi_\mathrm{rand} = 2 v \pi.
 * \f]
 *
 * Here, \f$u\f$ and \f$v\f$ denote independent uniform random numbers in the
 * range \f$\left[ 0, 1 \right]\f$. After that, a uniform random number
 * \f$W_\mathrm{rand}\f$ in the range \f$\left[ 0, W_\mathrm{max} \right]\f$ is
 * sampled. If
 *
 * \f[
 *      W_\mathrm{rand} \leq W \left( \theta_\mathrm{rand},
 * \varphi_\mathrm{rand} \right), \f]
 *
 * then the vector \f$\left( \theta_\mathrm{rand}, \varphi_\mathrm{rand}
 * \right)\f$ is accepted. If not, the vector is rejected and the algorithm
 * start from the beginning by sampling a new point on the sphere surface. Note
 * that the last inequality ensures that the probability of sampling a vector
 * \f$\left( \theta, \varphi \right)\f$ is proportional to \f$W \left( \theta,
 * \varphi \right)\f$. This algorithm obviously works best if the probability of
 * a vector being rejected is small, i.e. if \f$W_\mathrm{max}\f$ is a good
 * approximation of the real maximum of \f$W\f$, and if the distribution is not
 * too different from a uniform distribution \f$W = \left( 4 \pi
 * \right)^{-1}\f$. Note that rejection sampling does not require the
 * distribution to be normalized.
 *
 * In principle, rejection sampling proceeds in an infinite loop until a valid
 * vector has been found. Here, the loop is limited to a maximum number of trial
 * vectors \f$N_\mathrm{max}\f$. If this maximum number is reached, the Euler
 * angles \f$\left( 0, 0, 0 \right)\f$ are returned, which correspond to
 * spherical coordinates \f$\theta = 0\f$ and \f$\varphi = \pi/2\f$.
 */
class SphereRejectionSampler : public ReferenceFrameSampler {

public:
  /**
   * \brief Constructor
   *
   * \param dis \f$W \left( \theta, \varphi \right)\f$, probability distribution
   * in spherical coordinates. \param dis_max \f$W_\mathrm{max}\f$, upper limit
   * for the maximum of the probability distribution. \param seed Random number
   * seed. \param max_tri Maximum number of sampled points \f$\left(
   * \theta_\mathrm{rand}, \varphi_\mathrm{rand} \right)\f$ before the algorithm
   * terminates without success and returns \f$\left( 0, 0 \right)\f$ (default:
   * 1000).
   */
  SphereRejectionSampler(function<double(const double, const double)> dis,
                         const double dis_max, const unsigned long seed,
                         const unsigned int max_tri = 1000);

  /**
   * \brief Sample a random vector from probability distribution and record the
   * number of tries.
   *
   * \return std::pair which contains \f$N\f$, the number of tries that were
   * needed to find a valid vector, and the accepted vector \f$\left(
   * \theta_\mathrm{rand}, \varphi_\mathrm{rand}\right)\f$. Returns a std::pair
   * of \f$N_\mathrm{max}\f$ and \f$\left(0, 0 \right)\f$, if the maximum number
   * of trials \f$N_\mathrm{max}\f$ is reached by the algorithm and no random
   * vector was accepted.
   */
  pair<unsigned int, EulerAngles> sample() override;

protected:
  /**
   * \brief Sample polar angle of a uniformly randomly distributed point on a
   * sphere surface.
   *
   * \return \f$\theta_\mathrm{rand}\f$, random polar angle.
   */
  double sample_theta();

  /**
   * \brief Sample azimuthal angle of a uniformly randomly distributed point on
   * a sphere surface.
   *
   * \return \f$\varphi_\mathrm{rand}\f$, random azimuthal angle.
   */
  double sample_phi();

  /**
   * \brief Sample uniformly randomly distributed point on a sphere surface.
   *
   * \return \f$\left( \theta_\mathrm{rand}, \varphi_\mathrm{rand} \right)\f$,
   * random point on sphere surface.
   */
  CoordDir sample_theta_phi();

  function<double(const double, const double)>
      distribution; /**< \f$W \left( \theta, \varphi \right)\f$, (unnormalized)
                       probability distribution. */
  double distribution_maximum; /**< \f$W_\mathrm{max}\f$, maximum of
                                        probability distribution. */
  unsigned int max_tries;      /**< \f$N_\mathrm{max}\f$, maximum number of
                                        tries to find a random vector. */

  mt19937 random_engine; /**< Deterministic random number engine. */
  uniform_real_distribution<double>
      uniform_random; /**< Uniform distribution from which all random numbers
                         are derived here. */
};

} // namespace alpaca
