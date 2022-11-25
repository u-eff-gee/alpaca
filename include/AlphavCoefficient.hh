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

#include "FCoefficient.hh"
#include "KappaCoefficient.hh"
#include "StringRepresentable.hh"

/**
 * \brief Class for an \f$\alpha_\nu\f$ coefficient.
 * 
 * The coefficients \f$\alpha_\nu\f$ are a generalization of the \f$A_\nu\f$ coefficients that are 
 * encountered in the definition of the dir-dir correlation {Eq. (I-2) in \cite FaggHanna1959}.
 * To take into account information about the polarization of particles, the F coefficients which
 * constitute the \f$A_\nu\f$ coefficients are multiplied by their respective polarization 
 * coefficients \f$\kappa_\nu\f$ {Eq. (I-9) in \cite FaggHanna1959}:
 * 
 * \f[
 * 		\alpha_\nu \left( n \right)\left( L_n, L_n^\prime, j_n, j, \delta_n \right) = 
 * \f]
 * \f[
 * - \kappa_nu \left( L_n, L_n \right) F_\nu \left( L_n, L_n, j_n, j \right) 
 * \f]
 * \f[
 * + 2 \delta_n \kappa_\nu \left( L_n, L_n^\prime \right) F_\nu \left( L_n, L_n^\prime, j_n, j \right)
 * \f]
 * \f[
 * + \delta_n^2 \kappa_\nu \left( L_n^\prime, L_n^\prime \right) F_\nu \left( L_n^\prime, L_n^\prime, j_n, j \right).
 * \f]
 * 
 * See also the definition of the AvCoefficient class for more information.
 */

class AlphavCoefficient : public StringRepresentable{
public:

	/**
	 * \brief Constructor
	 * 
	 * \param two_nu \f$2 \nu\f$
	 * \param two_L Primary multipolarity \f$2 L\f$
	 * \param two_Lp Secondary multipolarity \f$2 L^\prime\f$
	 * \param two_jn Angular momentum quantum number \f$2 j_n\f$ of the initial or final state
	 * 	of a transition 
	 * \param two_j Angular momentum quantum number \f$2 j\f$ of the intermediate state
	 * 	of a transition 
	 */
	AlphavCoefficient(const int two_nu, const int two_L, const int two_Lp, const int two_jn, const int two_j);	

	/**
	 * \brief Return value of a specific \f$\alpha_\nu\f$ coefficient.
	 *
	 * \param delta Multipole mixing ratio \f$\delta\f$
	 *
	 * \return \f$\alpha_\nu \left( L, L^\prime, j_n, j, \delta_n \right)\f$
	 */
	double operator()(const double delta) const;

	string string_representation(const unsigned int n_digits = 0, const vector<string> variable_names = {}) const;

protected:
	const int two_nu;
	const int two_L;
	const int two_Lp;
	const int two_jn;
	const int two_j;

	const FCoefficient constant_f_coefficient, linear_f_coefficient, quadratic_f_coefficient;
	const KappaCoefficient constant_kappa_coefficient, linear_kappa_coefficient, quadratic_kappa_coefficient;

	double constant_coefficient;
	double linear_coefficient;
	double quadratic_coefficient;

};
