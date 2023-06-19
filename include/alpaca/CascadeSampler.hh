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

#include <memory>
#include <vector>

#include <gsl/gsl_math.h>

#include "alpaca/AngCorrRejectionSampler.hh"
#include "alpaca/AngularCorrelation.hh"
#include "alpaca/ReferenceFrameSampler.hh"

using std::shared_ptr;
using std::vector;

namespace alpaca {

/**
 * \brief Sample directions of emission from an arbitrarily long cascade of
 * angular correlations.
 *
 * This class can be used to sample a cascade of gamma rays in which the
 * direction of emission of a gamma ray depends on the direction of emission of
 * the previous gamma ray.
 *
 * It makes use of the fact that an arbitrarily long cascade can be built from
 * two-step dir-dir or pol-dir correlations, as they are implemented in the
 * alpaca code (the simulation of pol-dir cascades is limited, see below). The
 * most obvious statement that this is possible can probably be found in the
 * introduction of a review article by Krane, Steffen, and Wheeler \cite
 * KraneSteffenWheeler1973. The authors base their formalism on the notion of an
 * "oriented state", which is an excited state whose magnetic substates are
 * populated unequally. The three possibilities of creating an oriented state
 * listed by the authors include the "observation of preceding radiation" along
 * a given axis \f$\vec{\mathbf{k}}_0\f$. In a multi-step cascade, this is
 * exactly what happens: a gamma ray from a previous emission defines the
 * reference frame \f$\vec{\mathbf{k}}_0\f$ for the subsequent two-step process.
 *
 * Since the sampling of random vectors from two-step cascades can be handled by
 * AngCorrRejectionSampler, the task of the present class is to rotate the
 * directional vectors from AngCorrRejectionSampler, which are defined with
 * respect to the \f$z\f$ axis, into the reference frame given by the direction
 * of emission sampled for the previous step. For the rotations, Euler angles in
 * the x convention (see also the EulerAngleRotation class) are used.
 *
 * In order to set up a cascade between \f$n > 3\f$ states with \f$n-1\f$
 * transitions
 *
 * \f[
 *      j_1 \left( {L_1} \atop {L_1^\prime} \right) j_2 \left( {L_2} \atop
 * {L_2^\prime} \right) j_3 ... j_{n-1} \left( {L_{n-1}} \atop {L_{n-1}^\prime}
 * \right) j_n, \f]
 *
 * a std::vector of \f$n-2\f$ AngularCorrelation objects is passed to the
 * constructor of CascadeSampler, with the first element denoting the
 * first two transitions,
 *
 * \f[
 *      j_1 \left( {L_1} \atop {L_1^\prime} \right) j_2 \left( {L_2} \atop
 * {L_2^\prime} \right) j_3, \f]
 *
 * the second element describing the second and third transition,
 *
 * \f[
 *      j_2 \left( {L_2} \atop {L_2^\prime} \right) j_3 \left( {L_3} \atop
 * {L_3^\prime} \right) j_4, \f]
 *
 * and so on.
 * Note that the transition from the initial state to the intermediate state of
 * the two-step correlation \f$m\f$ is the same as the transition from the
 * intermediate state to the final state for the two-step correlation \f$m-1\f$
 * for all intermediate steps \f$1 < m < n\f$. This means that \f$n-3\f$
 * transitions appear twice in the AngularCorrelation objects.
 *
 * The first angular correlation is special in the sense that a reference frame
 * must be defined. There are two possibilities:
 *
 *  1. The cascade originates from an unoriented source. In that case, the first
 * gamma ray is emitted in a direction sampled from a uniform random
 * distribution on the surface of a sphere.
 *  2. The first gamma ray represents one of the mechanisms described by Krane,
 * Steffen, and Wheeler \cite KraneSteffenWheeler1973 that orients the second
 * state. For example, in a nuclear resonance fluorescence experiment, the first
 * transition could correspond to the capture of a beam photon.
 *
 * In the first case, it is usually desirable to return the direction of the
 * first gamma ray, since it belongs to the cascade. In the second case, the
 * user can decide whether the first gamma ray should be included in the output.
 *
 * Please note that the alpaca code implements dir-dir and pol-dir correlations,
 * but no pol-pol correlations, which would be required to simulate the
 * evolution of the gamma-ray polarization as the cascade progresses.
 * Consequently, only a single polarization measurement (i.e. a single instance
 * of an AngularCorrelation with fully specified electromagnetic transitions)
 * should be included in the cascade to avoid unphysical results. Again, in the
 * example of a nuclear resonance fluorescence experiment, this could be a
 * polarized photon beam that excites a nucleus, leading to a cascade of the
 * type
 *
 * \f[
 *      j_1 \left( {\vec{L_1}} \atop {L_1^\prime} \right) j_2 \left( {L_2} \atop
 * {L_2^\prime} \right) j_3 ... j_{n-1} \left( {L_{n-1}} \atop {L_{n-1}^\prime}
 * \right) j_n. \f]
 *
 * Since only the orientation of the polarization plane of the first transition
 * can be controlled by the user (only in case 2, where the second state is
 * oriented) by passing a set of Euler angles to the constructor of
 * CascadeSampler, all other transitions can be treated as if they were
 * rotationally symmetric around the \f$z\f$ axis. This means that if the single
 * polarized transition that may show up in the cascade is not the second
 * transition, the polarization plane can simply be sampled from a uniform
 * random distribution because there is nothing known about its orientation.
 * This implies that the relation between a vector \f$\left( \theta, \varphi
 * \right)\f$ in spherical coordinates and a set of Euler angles \f$\left( \Phi,
 * \Theta, \Psi \right)\f$ which, applied to the Cartesian vector \f$\left( 0,
 * 0, 1 \right)\f$ (i.e the \f$z\f$ axis), generates this vector, is unique: The
 * first Euler rotation around the \f$z\f$ axis by the angle \f$\Phi\f$ has no
 * effect on a rotationally symmetric body, therefore it can be neglected.
 * Internally, it is set to zero.
 * The orientation of the angular-correlation pattern along the last emitted
 * gamma ray can be achieved using the Euler angles \f$\Theta\f$ and \f$\Psi\f$
 * only.
 *
 * The CascadeSampler class uses a single random-number seed and a
 * single maximum number of iterations for all of its random-number sampling
 * members (\f$n-2\f$ instances of AngCorrRejectionSampler, and a single
 * instance of SphereRejectionSampler).
 */
class CascadeSampler {

public:
  explicit CascadeSampler(vector<shared_ptr<ReferenceFrameSampler>> cascade);

  /**
   * \brief Sample random gamma-ray directions from the cascade.
   *
   * \return List of reference frames in Euler angles. The first reference frame
   * describes the direction of emission of the first (depends on the setting of
   * return_first_direction) gamma ray in the cascade, the second pair describes
   * the second gamma ray, and so on.
   */
  vector<EulerAngles> operator()();

  [[nodiscard]] inline size_t size() const {
    return angular_correlation_samplers.size();
  }

protected:
  vector<shared_ptr<ReferenceFrameSampler>>
      angular_correlation_samplers; /**< List of AngCorrRejectionSamplers which
                                       are initialized on construction with
                                       AngularCorrelation objects. */
};

} // namespace alpaca
