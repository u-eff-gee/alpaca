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

#include <string>

#include "alpaca/FCoefficient.hh"
#include "alpaca/KappaCoefficient.hh"
#include "alpaca/StringRepresentable.hh"

using std::string;

namespace alpaca {

/**
 * \brief Class for an \f$\alpha_\nu\f$ coefficient.
 *
 * The coefficients \f$\alpha_\nu\f$ are a generalization of the \f$A_\nu\f$
 * coefficients that are encountered in the definition of the dir-dir
 * correlation {Eq. (I-2) in \cite FaggHanna1959}. To take into account
 * information about the polarization of particles, the F coefficients which
 * constitute the \f$A_\nu\f$ coefficients are multiplied by their respective
 * polarization coefficients \f$\kappa_\nu\f$ {Eq. (I-9) in \cite
 * FaggHanna1959}:
 *
 * \f[
 * 		\alpha_\nu \left( n \right)\left( L_n, L_n^\prime, j_n, j,
 * \delta_n \right) = \f] \f[
 * - \kappa_nu \left( L_n, L_n \right) F_\nu \left( L_n, L_n, j_n, j \right)
 * \f]
 * \f[
 * + 2 \delta_n \kappa_\nu \left( L_n, L_n^\prime \right) F_\nu \left( L_n,
 * L_n^\prime, j_n, j \right) \f] \f[
 * + \delta_n^2 \kappa_\nu \left( L_n^\prime, L_n^\prime \right) F_\nu \left(
 * L_n^\prime, L_n^\prime, j_n, j \right). \f]
 *
 * See also the definition of the AvCoefficient class for more information.
 */

class AlphavCoefficient : public StringRepresentable {
public:
  /**
   * \brief Constructor
   *
   * \param two_nu \f$2 \nu\f$
   * \param two_L Primary multipolarity \f$2 L\f$
   * \param two_Lp Secondary multipolarity \f$2 L^\prime\f$
   * \param two_jn Angular momentum quantum number \f$2 j_n\f$ of the initial or
   * final state of a transition \param two_j Angular momentum quantum number
   * \f$2 j\f$ of the intermediate state of a transition
   */
  AlphavCoefficient(int a_two_nu, int a_two_L, int a_two_Lp, int a_two_jn,
                    int a_two_j)
      : two_nu(a_two_nu), two_L(a_two_L), two_Lp(a_two_Lp), two_jn(a_two_jn),
        two_j(a_two_j),
        constant_f_coefficient(a_two_nu, a_two_L, a_two_L, a_two_jn, a_two_j),
        linear_f_coefficient(a_two_nu, a_two_L, a_two_Lp, a_two_jn, a_two_j),
        quadratic_f_coefficient(a_two_nu, a_two_Lp, a_two_Lp, a_two_jn,
                                a_two_j),
        constant_kappa_coefficient(a_two_nu, a_two_L, a_two_L),
        linear_kappa_coefficient(a_two_nu, a_two_L, a_two_Lp),
        quadratic_kappa_coefficient(a_two_nu, a_two_Lp, a_two_Lp),
        constant_coefficient(-constant_kappa_coefficient.get_value() *
                             constant_f_coefficient.get_value()),
        linear_coefficient(2. * linear_kappa_coefficient.get_value() *
                           linear_f_coefficient.get_value()),
        quadratic_coefficient(quadratic_kappa_coefficient.get_value() *
                              quadratic_f_coefficient.get_value()) {}

  /**
   * \brief Return value of a specific \f$\alpha_\nu\f$ coefficient.
   *
   * \param delta Multipole mixing ratio \f$\delta\f$
   *
   * \return \f$\alpha_\nu \left( L, L^\prime, j_n, j, \delta_n \right)\f$
   */
  inline double operator()(const double delta) const {
    return constant_coefficient + delta * linear_coefficient +
           delta * delta * quadratic_coefficient;
  }

  string string_representation(const int n_digits = 0,
                               const vector<string> variable_names = {}) const;

protected:
  int two_nu;
  int two_L;
  int two_Lp;
  int two_jn;
  int two_j;

  FCoefficient constant_f_coefficient, linear_f_coefficient,
      quadratic_f_coefficient;
  KappaCoefficient constant_kappa_coefficient, linear_kappa_coefficient,
      quadratic_kappa_coefficient;

  double constant_coefficient;
  double linear_coefficient;
  double quadratic_coefficient;
};

} // namespace alpaca
