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

#include <vector>

using std::vector;

#include <utility>

using std::pair;

#include "State.hh"
#include "StringRepresentable.hh"
#include "Transition.hh"

/**
 * \brief Abstract class for angular correlations of two photons.
 *
 * Base class for angular correlations of two photons, which defines the general
 * API. A gamma-gamma angular correlation \f$W \left( \theta, \varphi
 * \right)\f$ is defined by an oriented initial state and a set of cascade
 * steps, which are pairs of a transition and a state to which this transition
 * leads. The angular correlation is a function of two variables, the polar
 * angle \f$\theta\f$ and the azimuthal angle \f$\varphi\f$ in spherical
 * coordinates.
 */
class W_gamma_gamma : public StringRepresentable {

public:
  /**
   * \brief Constructor
   *
   * \param ini_sta Oriented intial state.
   * \param cas_ste Steps of the cascade which follow the exciation of the
   * initial state. Each step is a pair of a transition and the state which is
   * populated by that transition.
   */
  W_gamma_gamma(const State &ini_sta,
                const vector<pair<Transition, State>> cas_ste)
      : initial_state(ini_sta), cascade_steps(cas_ste),
        n_cascade_steps(cas_ste.size()) {}

  /**
   * \brief Destructor
   *
   * Virtual destructor needed to ensure that destructors of derived classes
   * W_dir_dir and W_pol_dir are called whenever a W_gamma_gamma goes out of
   * scope. Since AngularCorrelation, the main user interface, is always in
   * posession of a unique_ptr to one of the derived classes, not having this
   * destructor lead to memory leaks in the past.
   */
  virtual ~W_gamma_gamma(){};

  /**
   * \brief Call operator of the gamma-gamma angular correlation
   *
   * Returns the value of the angular correlation at a polar angle \f$\theta\f$
   * and an azimuthal angle \f$\varphi\f$ in spherical coordinates.
   *
   * \param theta Polar angle in spherical coordinates in radians
   * (\f$\theta \in \left[ 0, \pi \right]\f$).
   * \param phi Azimuthal angle in spherical coordinates in radians
   * (\f$\varphi \in \left[ 0, 2 \pi \right]\f$).
   *
   * \return \f$W_{\gamma \gamma} \left( \theta, \varphi \right)\f$
   */
  virtual double operator()(const double theta, const double phi) const = 0;

  /**
   * \brief Return an upper limit for possible values of the gamma-gamma angular
   * correlation.
   *
   * Some applications, for example the rejection-sampling (or 'accept-reject')
   * algorithm (see, e.g. Sec. 2.3 in \cite RobertCasella1999), which can be
   * used to sample random directions that are distributed according to a given
   * angular correlation, require an expression, or at least an estimate, for
   * the maximum absolute value of \f$W \left( \theta, \varphi \right)\f$, i.e.:
   *
   * \f[
   * 		\mathrm{max}_{\theta \in \left[ 0, \pi \right], \varphi \in
   * \left[ 0, 2\pi \right]} | W \left( \theta, \varphi \right) |. \f]
   *
   * If a useful upper limit estimate exists for a given angular correlation,
   * this function will return it. If no useful upper limit exists, or the
   * absolute value of \f$W\f$ does not have a limit, this function returns a
   * negative number.
   *
   * \return \f$\mathrm{max}_{\theta \in \left[ 0, \pi \right], \varphi \in
   * \left[ 0, 2\pi \right]} | W \left( \theta, \varphi \right) | \f$, or an
   * upper limit for this quantity. If no useful upper limit can be given or if
   * there is no limit, a negative number is returned.
   */
  virtual double get_upper_limit() const = 0;

  /**
   * \brief Return the initial state of the angular correlation.
   *
   * \return Initial state.
   */
  State get_initial_state() const { return initial_state; }

  /**
   * \brief Return the cascade steps.
   *
   * \return vector of Transition-State pairs.
   */
  vector<pair<Transition, State>> get_cascade_steps() const {
    return cascade_steps;
  }

  /**
   * @brief Calculate the normalization factor for the angular correlation.
   *
   * In this library, all angular correlations are supposed to be normalized to
   * \f$4\pi\f$ with respect to an integration over the entire surface of a
   * sphere, i.e.:
   *
   * \f[
   *  \int_0^{2 \pi} \mathrm{d} \varphi \int_0^\pi W \left( \theta, \varphi
   * \right) \sin \left( \theta \right) \mathrm{d} \theta = 4\pi.
   * \f]
   *
   * Considering that the integral over the entire range of \f$\theta\f$ is zero
   * except for the zero-order Legendre polynomial
   *
   * \f[
   *  \int_0^\pi P_n \left[ \cos \left( \theta \right) \right] \mathrm{d} \theta
   * = \begin{cases} 1,~~~n=0 \\ 0,~~~\mathrm{else} \end{cases}, \f]
   *
   * and that the integral over the entire range of \f$\varphi\f$ is zero for
   * the angular dependence of the polarization-dependent term of the pol-dir
   * correlation
   *
   * \f[
   *  \int_0^{2\pi} \cos \left( 2 \varphi \right) = 0,
   * \f]
   *
   * it is clear that only the zero-order terms of the dir-dir correlation need
   * to be considered to determine the normalization factor (the value of \f$4
   * \pi\f$ comes from the integration over the azimuthal angle when the
   * integration over the polar angle does not vanish).
   *
   * For the most general angular correlation with unobserved intermediate
   * \f$\gamma\f$ rays, the integral over the unnormalized distribution {i.e.
   * using the canonical \f$A_\nu\f$, \f$\alpha_\nu\f$, \f$E_\nu\f$, and
   * \f$U_\nu\f$ coefficients as in Refs. \cite FaggHanna1959 [Eqs. I-1, I-2,
   * I-7, I-8, and I-9] or \cite Iliadis2021 [Eqs. (3-7) and (21-24)]} results
   * in [for explicit examples, see Sec. 4 in Ref. \cite Iliadis2021 or Ref.
   * \cite AjzenbergSelove1960]:
   *
   * \f[
   *  \int_0^{2 \pi} \mathrm{d} \varphi \int_0^\pi W_\mathrm{unnormalized}
   * \left( \theta, \varphi \right) \sin \left( \theta \right) \mathrm{d} \theta
   * = \f]
   *
   * \f[
   * \left[ F_0 \left(
   *  L_1, L_1, j_1, J \right) + F_0 \left( L_1, L_1^\prime, j_1, J \right)
   *  \delta_1 + F_0 \left( L_1^\prime, L_1^\prime, j_1, J \right) \delta_1^2
   *  \right]
   * \f]
   *
   * \f[
   *  \times \left[ 1 + \delta_2^2
   *  \right] \times ... \times \left[ 1 + \delta_{n-2}^2
   *  \right]
   * \f]
   *
   * \f[
   *  \times \left[ F_0 \left(
   *  L_{n-1}, L_{n-1}, j_n, J \right) + F_0 \left( L_{n-1}, L_{n-1}^\prime,
   * j_n, J \right) \delta_{n-1} + F_0 \left( L_{n-1}^\prime, L_{n-1}^\prime,
   * j_n, J \right) \delta_{n-1}^2 \right]
   * \f].
   *
   * Please note that the unobserved \f$\gamma\f$-ray transitions mix
   * incoherently with the others, so the mixing ratios of those can be taken
   * into account in a more straightforward way than for correlated \f$\gamma\f$
   * rays [see, e.g., Eq. (21) in Ref. \cite Iliadis2021].
   *
   * The expression above is often simplified by using special properties of
   * zero-order F coefficients. Firstly, the constant and quadratic terms in the
   * mixing ratio all evaluate to 1 [see, e.g., Eq. (9) in Ref. \cite
   * FerentzRosenzweig1955]:
   *
   * \f[
   *  F_0 \left( L, L, j, J \right) = 1
   * \f]
   *
   * if they describe a valid transition (see also below).
   * Secondly, all linear terms vanish [See, e.g., Eq. (10) in Ref. \cite
   * Iliadis2021. The first part of the inequality cannot be fulfilled if \f$L
   * \neq L^\prime\f$ and \f$n=0\f$ in the notation of that article.]
   *
   * \f[
   *  F_0 \left( L, L^\prime, j, J \right) = 0~~\forall~~L \neq L^\prime.
   * \f]
   *
   * Therefore, the expression above becomes:
   *
   * \f[
   *  \int_0^{2 \pi} \mathrm{d} \varphi \int_0^\pi W_\mathrm{unnormalized}
   * \left( \theta, \varphi \right) \sin \left( \theta \right) \mathrm{d} \theta
   * = \f]
   *
   * \f[
   * \left( 1 + \delta_1^2 \right) \left( 1 + \delta_2^2 \right) \times ...
   * \times \left( 1 + \delta_{n-2}^2 \right) \left( 1 +
   * \delta_{n-1}^2 \right). \f]
   *
   * Simplified expressions like the one above can be found in recent articles,
   * for example Refs. \cite Iliadis2021 and \cite Zilges2022.
   * However, the alpaca library uses the original expression with the F
   * coefficients for the following reason: Imagine a cascade that includes a
   * transition where only a single multipolarity is possible (for example a
   * transition involving a spin-zero state). In this case, the quadratic
   * zero-order F coefficients would evaluate to zero as well since it violates
   * the triangle inequality. Therefore, the normalization factor would be
   * independent of the mixing ratio for that transition. Knowing this, Eq. (12)
   * in \cite Iliadis2021 could be misleading in the sense that it implies that
   * the numerical result of the canonical angular-correlation expression always
   * need to be divided by the mixing ratios. If a user forgets that a
   * transition only admits a single multipolarity and sets the mixing ratio to
   * a nonzero value, this would break the normalization.
   * I chose to mention this in detail here, since alpaca versions before
   * version 1.0.5 implemented the normalization factor in a way that the
   * normalization could fail.
   *
   * @return Normalization factor.
   */
  // virtual double calculate_normalization_factor() const = 0;

  /**
   * @brief Return string representation of gamma-gamma angular correlation.
   *
   * \param n_digits Determines whether the expression
   * should be evaluated numerically, (n_digits > 0). If yes, indicates how many
   * digits should be displayed (default: 0).
   * \param variable_names Names for the variables of a function (default: {}
   * i.e. use default names).
   *
   * \return String representation.
   */
  virtual string string_representation(
      const unsigned int n_digits = 0,
      const vector<string> variable_names = {}) const override = 0;

protected:
  State initial_state; /**< Initial state */
                       /**
                        * Steps of the gamma-ray cascade following an excitation.
                        * Each step consists of an electromagnetic transition and a state which is
                        * populated by                      that transition.
                        */
  vector<pair<Transition, State>> cascade_steps;

  double normalization_factor;  /**< Normalization factor for the angular
                                   distribution */
  const size_t n_cascade_steps; /**< Number of transitions in the cascade. */
  int two_nu_max; /**< Maximum value of \f$2 \nu\f$ for which the coefficients
                     do not vanish */
  int nu_max; /**< Maximum value of \f$\nu\f$ for which the coefficients do not
                 vanish */
};