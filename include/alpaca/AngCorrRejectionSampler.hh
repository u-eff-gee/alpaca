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

#include "alpaca/AngularCorrelation.hh"
#include "alpaca/SphereRejectionSampler.hh"

namespace alpaca {

/**
 * \brief Sample directions in spherical coordinates from an angular
 * correlation.
 *
 * Although SphereRejectionSampler already accepts an arbitrary function of
 * \f$\theta\f$ and \f$\varphi\f$, a function object of the class
 * AngularCorrelation cannot be passed this way. Therefore, the present class
 * was derived. For more information, see the base class.
 */
class AngCorrRejectionSampler : public SphereRejectionSampler {

public:
  /**
   * \brief Constructor
   *
   * In contrast to the base class, the upper limit \f$W_\mathrm{max}\f$ does
   * not have to be provided explicitly. The member function
   * AngularCorrelation::get_upper_limit() is called instead.
   *
   * \param w \f$W \left( \theta, \varphi \right)\f$, angular correlation
   * \param seed Random number seed.
   * \param max_tri Maximum number of sampled points
   * \f$\left( \theta_\mathrm{rand}, \varphi_\mathrm{rand} \right)\f$
   * before the algorithm terminates without success and returns \f$\left( 0, 0
   * \right)\f$.
   */
  AngCorrRejectionSampler(AngularCorrelation &w, const int seed,
                          const unsigned int max_tri = 1000);
};

} // namespace alpaca
