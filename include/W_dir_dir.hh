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

    Copyright (C) 2021, 2022 Udo Friman-Gayer
*/

#pragma once

#include "AvCoefficient.hh"
#include "UvCoefficient.hh"
#include "W_gamma_gamma.hh"

/**
 * \brief Class for a direction-direction (dir-dir) correlation
 *
 * Evaluates the angular correlation between two photons which are emitted from an oriented nuclear
 * state.
 * Without loss of generality, the orientation is assumed to be in the z direction.
 * It may result from a nuclear reaction of the type 
 * \f$\left( \gamma, \gamma^\prime \right)\f$, in which a nucleus in an initial state with the 
 * total angular momentum quantum number ('spin') \f$j_1\f$ absorbs a photon from a beam 
 * that travels in positive z direction to render it in an 
 * excited state with the spin \f$j_2\f$, 
 * and emits another photon in a transition to a lower-lying state with the spin \f$j_n\f$.
 * An arbitrary number of intermediate states \f$j_i\f$ with labels \f$1 < i < n\f$ may be
 * populated in the decay process.
 * It is assumed that only the two photons associated with the transitions 
 * \f$j_1 \to j_2\f$ and \f$j_{n-1} \to j_n\f$ are observed.
 * Furthermore, the restriction to a dir-dir correlation implies that only information about the
 * direction of the photons is available, but not about their polarization.
 * In particular, this means that the dir-dir correlation does not depend on the parity
 * quantum number of the involved nuclear states, which translates into an independence of the
 * azimuthal angle of the emitted photon.
 * 
 * The class uses Eqs. (I-1), (I-1') and (I-2) of Ref. \cite FaggHanna1959 with a similar notation.
 * Note that Eq. (I-2) already includes the assumption that only two multipoles contribute
 * to any of the two transitions, and that the spin of the intermediate states is known.
 * The process can be denoted as {this notation is similar to the one used by Biedenharn 
 * \cite AjzenbergSelove1960}:
 * 
 * \f[
 * 		j_1 \left( \begin{array}{c} L_1 \\ L_1^\prime \end{array} \right) j_2 \left( \begin{array}{c} L_2 \\ L_2^\prime \end{array} \right) ... j_n,
 * \f]
 * 
 * where the dots may represent an arbitrary number of intermediate transitions.
 * The entire sequence of transitions, whose first and last transition are observed, is called
 * a cascade.
 * 
 * The dir-dir correlation is normalized to \f$4 \pi\f$ here compared to Ref. \cite FaggHanna1959
 * {see below Eq. (I-2) therein} by dividing through the \f$\left( 1 + \delta^2 \right)\f$ factors.
 */
class W_dir_dir : public W_gamma_gamma{
public:
	/**
	 * \brief Constructor
	 * 
	 * \param ini_sta Oriented intial state.
	 * \param cas_ste Steps of the cascade which follow the exciation of the initial state.
	 * Each step is a pair of a transition and the state which is populated by that transition.
	 */
	W_dir_dir(const State &ini_sta, const vector<pair<Transition, State>> cas_ste);

	/**
	 * \brief Return value of the dir-dir correlation at an angle \f$\theta\f$
	 * 
	 * For a two-step cascade, it is given by the following expression 
	 * {Eqs. (I-1) and (I-2) in Ref. \cite FaggHanna1959}:
	 * 
	 * \f[
	 * 		W \left( \theta, \varphi \right) = W \left( \theta \right) = \sum_\nu A_v \left( L_1, L_1^\prime, j_1, j_2 \right) A_v \left( L_2, L_2^\prime, j_3, j_2 \right) P_\nu \left[ \cos \left( \theta \right) \right].
	 * \f]
	 * 
	 * The first transition proceeds via the multipolarities \f$L_1\f$ and \f$L_1^\prime\f$,
	 * while the second has the multipolarities \f$L_2\f$ and \f$L_2^\prime\f$.
	 * The symbol \f$P_v\f$ denotes a Legendre polynomial of the degree \f$\nu\f$.
	 * It incorporates the dependence on the polar angle \f$\theta\f$, while the dependence on the
	 * nuclear properties is encoded in the \f$A_\nu\f$ coefficients.
	 * 
	 * For a cascade with more than 2 steps, the so-called \f$U_\nu\f$ coefficients enter the 
	 * previous equation to take into account the deorientation due to the addition decays
	 * {Eq. (I-1') in Ref. \cite FaggHanna1959}:
	 * 
	 * \f[
	 * 		W \left( \theta \right) = \sum_\nu A_v \left( L_1, L_1^\prime, j_1, j_2 \right) U_\nu \left( j_2, L_2, L_2^\prime, j_3 \right) ... A_v \left( L_{n-1}, L_{n-1}^\prime, j_n, j_{n-1} \right) P_\nu \left[ \cos \left( \theta \right) \right].
	 * \f]
	 * 
	 * \param theta Polar angle between the direction of the incoming and 
	 * the outgoing photon in radians.
	 * 
	 * \return \f$W \left( \theta \right)\f$
	 */
	double operator()(const double theta) const;

	/**
	 * \brief Return value of the dir-dir correlation at an angle \f$\theta\f$
	 * 
	 * This call operator accepts the angle \f$\varphi\f$ as a second argument, although the 
	 * direction-direction correlation is independent of the azimuthal angle.
	 * Internally, the method calls the single-argument call operator.
	 * 
	 * \param theta Polar angle between the direction of the incoming and 
	 * the outgoing photon in radians.
	 * \param phi Azimuthal angle between the polarization axis of the first photon and the direction of the
	 * outgoing photon in radians. Note that this argument has no influence on the result.
	 * 
	 * \return \f$W \left( \theta \right)\f$
	 */
	double operator()(const double theta, const double) const override {
		return operator()(theta);
	}
	
	/**
	 * \brief Return upper limit for the dir-dir correlation.
	 * 
	 * For the dir-dir correlation, the maximum possible value cannot be calculated analytically,
	 * but an upper limit can be given.
	 * In the following, the angular correlation without unobserved intermediate transitions is
	 * considered without loss of generality.
	 * Start from the exact expression for the maximum absolute value of \f$W\f$ 
	 * {Eq. (I-1') in Ref. \cite FaggHanna1959} and apply the triangle inequality:
	 * 
	 * \f[
	 * 		\mathrm{max}_{\theta \in \left[ 0, \pi \right], \varphi \in \left[ 0, 2\pi \right]} 
	 * | W \left( \theta, \varphi \right) | = 
	 * \f]
	 * \f[
	 * 		= \mathrm{max}_{\theta \in \left[ 0, \pi \right], \varphi \in \left[ 0, 2\pi \right]} 
	 * | \sum_\nu A_v \left( L_1, L_1^\prime, j_1, j_2 \right) A_v \left( L_2, L_2^\prime, j_3, j_2 \right) P_\nu \left[ \cos \left( \theta \right) \right] |
	 * \f]
	 * \f[
	 * 		\leq 
	 * \sum_\nu | A_v \left( L_1, L_1^\prime, j_1, j_2 \right) A_v \left( L_2, L_2^\prime, j_3, j_2 \right) | \mathrm{max}_{\theta \in \left[ 0, \pi \right]} | P_\nu \left[ \cos \left( \theta \right) \right] |.
	 * \f]
	 * 
	 * The last line summarizes the result of three steps,
	 * 
	 * \f[
	 * 		\mathrm{max}_x | \sum_i a_i f_i \left( x \right) | 
	 * 		\leq \mathrm{max}_x \sum_i | a_i f_i \left( x \right) | 
	 * 		= \mathrm{max}_x \sum_i | a_i | | f_i \left( x \right) | 
	 * 		\leq \sum_i | a_i | | \mathrm{max}_x f_i \left( x \right) |,
	 * \f]
	 * 
	 * where the $a_i$ are real coefficients and the \f$f_i\f$ real-valued function of a 
	 * scalar variable \f$x\f$.
	 * 
	 * For the Legendre polynomials \f$P_\nu\f$, the following inequality holds [See, e.g. 
	 * Eq. (18.14.1) in \cite DLMF2020, which gives an upper limit for the more general Jacobi
	 * polynomials.\
	 * The relation between Jacobi and Legendre polynomials is given in Eq. (18.7.9) therein]:
	 * 
	 * \f[
	 * 		|P_\nu \left[ \cos \left( \theta \right) \right]| \leq P_\nu \left( 1 \right) = 1.
	 * \f]
	 * 
	 * Therefore, this function returns
	 * 
	 * \f[
	 * 		\sum_\nu | A_v \left( L_1, L_1^\prime, j_1, j_2 \right) A_v \left( L_2, L_2^\prime, j_3, j_2 \right) | \geq \mathrm{max}_{\theta \in \left[ 0, \pi \right], \varphi \in \left[ 0, 2\pi \right]} 
	 * | W \left( \theta, \varphi \right) |
	 * \f]
	 * 
	 * as an upper limit.
	 * 
	 * \return \f$\mathrm{max}_{\theta \in \left[ 0, \pi \right], \varphi \in \left[ 0, 2\pi \right]} 
	 * | W \left( \theta, \varphi \right) | \f$, or an upper limit for this quantity.
	 * If no useful upper limit can be given or if there is no limit, a negative number is returned.
	 */
	double get_upper_limit() const override;
	
	/**
	 * \brief Return \f$\nu_\mathrm{max}\f$
	 */
	int get_nu_max() const { return nu_max; };

	/**
	 * \brief Return \f$2 \nu_\mathrm{max}\f$
	 */
	int get_two_nu_max() const { return two_nu_max; };

	double get_normalization_factor() const { return normalization_factor; };

	/**
	 * \brief Return \f$U_\nu\f$ coefficients for the dir-dir correlation.
	 * 
	 * This function returns an \f$\nu_\mathrm{max} / 2 \times \left(n-2\right) f$ array of 
	 * \f$U_\nu\f$ coefficient (UvCoefficient) objects as opposed to 
	 * W_dir_dir::get_Uv_coefficient_products() and 
	 * W_dir_dir::calculate_expansion_coefficients_Uv() (protected method), which return lists of 
	 * the products of the columns as a \f$\nu_\mathrm{max} / 2 \times 1\f$ array.
	 * 
	 * \return \f$U_\nu\f$ coefficients, sorted by \f$\nu_\mathrm{max}\f$ (first index) and 
	 * the cascade step number (second index, runs from \f$2\f$ to \f$n-1\f$).
	 */
	vector<vector<UvCoefficient>> get_Uv_coefficients() const { return uv_coefficients; };

	/**
	 * \brief Return products of \f$U_\nu\f$ coefficients to be inserted between the \f$A_\nu\f$/\f$/alpha_\nu\f$ coefficients .
	 * 
	 * As opposed to W_dir_dir::get_Uv_coefficients() (see also the documentation of that method),
	 * this function returns products of \f$U_\nu\f$ coefficients at a given value of \f$\nu\f$ 
	 * for all intermediate cascade steps.
	 * These are the values that need to be inserted between the \f$A_\nu\f$ or \f$\alpha_\nu\f$
	 * coefficients in an angular correlation to take into account the unobserved intermediate 
	 * steps [see, e.g. Eqs. (22) and (23) in Ref. \cite Iliadis2021].
	 * As opposed to W_dir_dir::calculate_expansion_coefficients_Uv(), this function does not 
	 * recalculate anything and it is meant for public access.
	 * 
	 * \return List of products of \f$U_\nu\f$ coefficients for all values of \f$\nu\f$.
	 */
	vector<double> get_Uv_coefficient_products() const { return uv_coefficient_products; };

	string string_representation(const unsigned int n_digits = 0, const vector<string> variable_names = {}) const override;

protected:
	/**
	 * \brief Calculate products of \f$U_\nu\f$ coefficients for the dir-dir correlation.
	 * 
	 * See also the documentation of W_dir_dir::calculate_expansion_coefficients().
	 * 
	 * \return \f$U_\nu (2) ... U_\nu (n-1)\f$ for \f$ \nu \in \lbrace 0, ..., \nu_\mathrm{max} \rbrace,~ \nu~\mathrm{even} \f$ sorted by increasing values
	 *  of \f$\nu\f$ in a std::vector.
	 */
	vector<double> calculate_expansion_coefficients_Uv();

	/**
	 * \brief Get the maximum value \f$\nu_\mathrm{max}\f$ for which the product of coefficients
	 * \f$A_\nu\f$ and \f$U_\nu\f$ is nonzero.
	 * 
	 * This limits the sum over \f$\nu\f$ in the definition of the direction-direction correlation.
	 * The respective upper limits are given in the documentation of the AvCoefficient and 
	 * UvCoefficient classes.
	 * 
	 * This function calls the function W_dir_dir::calculate_two_nu_max_Av(), 
	 * which calculates the upper limits for the \f$A_\nu\f$ coefficients, and, if there are 
	 * unobserved intermediate transitions, the function W_dir_dir::calculate_two_nu_max_Uv() 
	 * to do the same for the \f$U_\nu\f$ coefficients.
	 * The net limit for \f$\nu\f$ will be the more restrictive of both.
	 * 
	 * Note that the same limit also applies for polarization-direction correlations.
	 * Although they make use of more general coefficients than the \f$A_\nu\f$, they still depend
	 * on the same F coefficients {compare, e.g., Eq. (I-2) and (I-9) in \cite FaggHanna1959}.
	 *
	 * \param ini_to_int Transition from initial state to intermediate state
	 * \param int_state Intermediate state
	 * \param int_to_fin Transition from intermediate state to final state
	 * 
	 * \return \f$2\nu_\mathrm{max}\f$
	 */
	int calculate_two_nu_max() const;

	/**
	 * \brief Get the maximum value \f$\nu_\mathrm{max}\f$ for which the product of coefficients
	 * \f$A_\nu\f$ is nonzero.
	 * 
	 * See also the documentation of W_dir_dir::calculate_two_nu_max().
	 * 
	 * \return \f$2\nu_\mathrm{max}\f$, restriction due to properties of the \f$A_\nu\f$ coefficients
	 */
	int calculate_two_nu_max_Av() const;

	/**
	 * \brief Get the maximum value \f$\nu_\mathrm{max}\f$ for which the product of coefficients
	 * \f$U_\nu\f$ is nonzero.
	 * 
	 * See also the documentation of W_dir_dir::calculate_two_nu_max().
	 * 
	 * \return \f$2\nu_\mathrm{max}\f$, restriction due to properties of the \f$U_\nu\f$ coefficients
	 */
	int calculate_two_nu_max_Uv() const;

	/**
	 * \brief Calculate the set of expansion coefficients for the dir-dir correlation.
	 * 
	 * The sum over \f$\nu\f$ in Eqs. (I-1) and (I-1') of \cite FaggHanna1959 contains products of 
	 * \f$A_\nu\f$ and potentially also \f$U_\nu\f$ coefficients.
	 * 
	 * This function calls the function W_dir_dir::calculate_expansion_coefficients_Av() to 
	 * calculate the products of \f$A_\nu\f$ coefficients from \f$\nu = 0\f$ up to, and including, 
	 * a maximum value of \f$\nu_\mathrm{max}\f$ for a given set of quantum numbers.
	 * If the cascade contains more than two transitions, the corresponding products of 
	 * \f$U_\nu\f$ coefficients for the same values of \f$\nu\f$ are calculated by the 
	 * W_dir_dir::calculate_expansion_coefficients_Av() function.
	 * This function merges the output of the two functions.
	 * 
	 * \return \f$A_\nu (1) A_\nu (2)\f$ (\f$n=3\f$) or \f$A_\nu (1) U_\nu (2) ... A_\nu (n)\f$ (\f$n>3\f$) for \f$ \nu \in \lbrace 0, ..., \nu_\mathrm{max} \rbrace,~ \nu~\mathrm{even} \f$ sorted by increasing values
	 *  of \f$\nu\f$ in a std::vector.
	 */
	vector<double> calculate_expansion_coefficients();
	
	/**
	 * \brief Calculate products of \f$A_\nu\f$ coefficients for the dir-dir correlation.
	 * 
	 * See also the documentation of W_dir_dir::calculate_expansion_coefficients().
	 * 
	 * \return \f$A_\nu (1) A_\nu (n)\f$ for \f$ \nu \in \lbrace 0, ..., \nu_\mathrm{max} \rbrace,~ \nu~\mathrm{even} \f$ sorted by increasing values
	 *  of \f$\nu\f$ in a std::vector.
	 */
	vector<double> calculate_expansion_coefficients_Av();

	/**
	 * \brief Calculate the normalization factor for the angular correlation.
	 * 
	 * As stated below Eqs. (I-2) and (I-1') in Ref. \cite FaggHanna1959, the \f$A_\nu\f$ 
	 * and \f$U_\nu\f$ coefficients are 
	 * normalized to \f$1 + \delta^2\f$, where \f$\delta\f$ is the multipole mixing ratio that
	 * they contain.
	 * 
	 * In order to normalize the angular correlations to \f$4 \pi\f$, it must be multiplied by:
	 * 
	 * \f[
	 * 		\prod_{i = 1}^{n-1} \frac{1}{1+\delta_i^2},
	 * \f]
	 * 
	 * where \f$\delta_i\f$ is the multipole mixing ratio of the \f$i\f$-th transition.
	 * 
	 * \return \f$\prod_{i} \left( 1+\delta_i^2 \right)^{-1}\f$
	 */
	double calculate_normalization_factor() const;

	vector<AvCoefficient> av_coefficients_excitation; /**< Vector of AvCoefficient objects for the excitation */
	vector<AvCoefficient> av_coefficients_decay; /**< Vector of AvCoefficient objects for the decays */
	vector<vector<UvCoefficient>> uv_coefficients; /**< Vector of the UvCoefficient objects */
	vector<double> uv_coefficient_products; /**< Vector of products of \f$U_\nu\f$ coefficients */
	vector<double> expansion_coefficients; /**< Vector to store expansion coefficients */

};
