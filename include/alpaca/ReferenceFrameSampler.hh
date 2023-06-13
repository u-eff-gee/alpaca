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

#include <utility>

#include "alpaca/EulerAngleRotation.hh"

using std::pair;

using alpaca::EulerAngles;

namespace alpaca {

/**
 * @brief Abstract class for sampling an arbitrarily oriented reference frame.
 *
 * It is assumed that derived classes sample a reference frame with the Euler
 * angles \f$\Phi\f$, \f$\Theta\f$ and \f$\Psi\f$ using an algortihm that may
 * require a number of tries before a valid vector is found
 * [e.g. rejection sampling (see, e.g. Sec. 2.3 in Ref. \cite RobertCasella1999
 * and the *RejectionSampler classes for more information].
 *
 * This class defines the abstract interface for the samplers, including a
 * method to report the number of tries until a valid vector was found
 *
 */
class ReferenceFrameSampler {
public:
  virtual ~ReferenceFrameSampler() = default;
  ReferenceFrameSampler() = default;
  ReferenceFrameSampler(const ReferenceFrameSampler &) = default;
  ReferenceFrameSampler &operator=(const ReferenceFrameSampler &) = default;
  ReferenceFrameSampler(ReferenceFrameSampler &&) = default;
  ReferenceFrameSampler &operator=(ReferenceFrameSampler &&) = default;

  /**
   * \brief Sample a random reference frame and record the
   * number of tries.
   * Abstract function to be overriden by derived classes.
   *
   * \return std::pair which contains \f$N\f$, the number of tries that were
   * needed to find a valid vector, and the accepted reference frame \f$\left(
   * \Phi_\mathrm{rand}, \Theta_\mathrm{rand}, \Psi_\mathrm{rand}\right)\f$.
   * Returns a std::pair of \f$N_\mathrm{max}\f$ and \f$\left(0, 0, 0
   * \right)\f$, if the maximum number of trials \f$N_\mathrm{max}\f$ is reached
   * by the algorithm and no random reference frame was accepted.
   */
  virtual pair<unsigned int, EulerAngles> sample() = 0;
  /**
   * \brief Sample a random reference frame.
   *
   * \return Accepted reference frame \f$\left( \theta_\mathrm{rand},
   * \varphi_\mathrm{rand}\right)\f$. Returns \f$\left( 0, 0 \right)\f$ if the
   * maximum number of trials \f$N_\mathrm{max}\f$ is reached by the algorithm
   * and no random vector was accepted.
   */
  EulerAngles operator()();

  /**
   * \brief Estimate the efficiency of sampling for the given
   * distribution.
   *
   * The efficiency \f$\epsilon\f$ is estimated by sampling \f$n\f$ reference
   * frames from the distribution and calculating the average number of trials
   * \f$\langle N \rangle\f$ for this set.
   *
   * \param n_tries \f$n\f$, number of sampled vectors.
   *
   * \return Estimate for \f$\epsilon\f$ from the \f$n\f$ samples.
   */
  double estimate_efficiency(const unsigned int n_tries);
};

} // namespace alpaca
