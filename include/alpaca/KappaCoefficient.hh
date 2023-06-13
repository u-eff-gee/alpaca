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

#include "alpaca/StringRepresentable.hh"

namespace alpaca {

/**
 * \brief Class for a polarization coefficient
 *
 * The coefficients \f$\kappa_\nu\f$ are introduced in Eq. (I-4) of Ref. \cite
 *FaggHanna1959 to take into account polarization information, which is probably
 *why they are called polarization coefficients in \cite Kneissl1996. They are
 *defined as {Eq. (I-7) in \cite FaggHanna1959}:
 *
 * \f[
 *      \kappa_\nu \left( L, L^\prime \right) = - \sqrt{ \frac{ \left( \nu - 2
 *\right)! }{ \left( \nu + 2 \right)! } } \frac{ \left( \begin{array}{ccc}
 *			    L & L^\prime & \nu \\
 *			    1 & 1 & 2
 *		    \end{array}
 *	    \right)
 * }{
 *      \left(
 *		    \begin{array}{ccc}
 *			    L & L^\prime & \nu \\
 *			    1 & -1 & 0
 *		    \end{array}
 *	    \right).
 * }
 * \f]
 *
 * The definition includes two Clebsch Gordan (CG) coefficients.
 * For the CG coefficient in the denominator, divisions by zero are anticipated
 *if the triangle inequality {see, e.g., Eq. I.2.3.ii in \cite Messiah19622}
 *
 * \f[
 *      \left| L - L^\prime \right| < \nu < L + L^\prime
 * \f]
 *
 * is not fulfilled.
 * This is only a computational issue, since the ratio of all possible CG
 *coefficients encountered in the present formalism can be shown to have no pole
 *within the domain of possible \f$L\f$, \f$L^\prime\f$ and \f$\nu\f$ {see,
 *e.g., Eq. (85a) in \cite AjzenbergSelove1960}: The \f$\kappa_\nu\f$ always
 *appear in as a factor in front of an \f$F_\nu\f$ coefficient in terms like
 *{see, e.g., Eq. (I-9) in \cite FaggHanna1959}:
 *
 * \f[
 *      \kappa_\nu \left( L L^\prime \right) F_\nu \left( L, L^\prime, ...
 *\right). \f]
 *
 * The corresponding \f$F_\nu\f$ coefficient contains the same CG coefficient as
 *in the denominator of the \f$\kappa_\nu\f$ coefficient, which cancels the pole
 *(see the FCoefficient class). It also contains a Wigner-6j symbol which
 *imposes the triangle inequality above again. Despite this issue, the code uses
 *the CG coefficients since they allow the coefficients \f$\kappa_\nu\f$ to be
 *written in a compact form. Due to the considerations above, divisions by zero
 *are avoided by checking the triangle inequality and returning \f$\kappa_\nu =
 *0\f$ if it is violated.
 *
 * One will also run into problems for \f$\nu = 0\f$,
 * since two different angular momenta \f$L\f$ and \f$L^\prime\f$ can not be
 *coupled to zero, and the behavior of the factorial for negative
 *numbers is undefined (The latter problem persists for \f$\nu = 1\f$. However,
 *this case is not encountered here). Reference \cite FaggHanna1959 and Ref.
 *\cite LLoyd1952, which is cited by the former for \f$\kappa_\nu\f$, only give
 *values for the coefficients with \f$\nu > 2\f$ [see Tab. II(b) \cite
 *FaggHanna1959 and Tab. III \cite LLoyd1952, respectively, and note the
 *different normalization]. Equation (I-8) of Ref. \cite FaggHanna1959 mentions
 *a restriction to even, nonzero values of \f$\nu\f$ when the dir-dir
 *correlation is extended to include polarization, but does not attribute this
 *to the \f$\kappa_\nu\f$ coefficients. In the present implementation, it was
 *decided to throw an exception when \f$\kappa_{\nu < 2}\f$ is requested.
 */
class KappaCoefficient : public StringRepresentable {
public:
  /**
   * \brief Constructor, value of a specific polarization coefficent
   *
   * \f$\kappa_\nu \left( L, L^\prime \right) \f$, in particular
   * \f$\kappa_\nu = 0\f$ if the triangle inequality between \f$L\f$,
   * \f$L^\prime\f$, and \f$\nu\f$ is not fulfilled.
   *
   * \param two_nu \f$2 \nu\f$
   * \param two_L \f$2 L\f$
   * \param two_Lp \f$2 L^\prime \f$
   *
   * \exception std::invalid_argument If \f$\nu < 2\f$, i.e. \f$2\nu < 4\f$
   */
  KappaCoefficient(const int two_nu, const int two_L, const int two_Lp);

  double get_value() const { return value; };

  string string_representation([[maybe_unused]] const int n_digits = 0,
                               vector<string> variable_names = {}) const;

protected:
  int two_nu;
  int two_L;
  int two_Lp;

  double value;
};

} // namespace alpaca
