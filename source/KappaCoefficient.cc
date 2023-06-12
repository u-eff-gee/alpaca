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

#include <cmath>

#include <stdexcept>

using std::invalid_argument;

#include <string>

using std::to_string;

#include <gsl/gsl_sf.h>

#include "alpaca/KappaCoefficient.hh"
#include "alpaca/TestUtilities.hh"

KappaCoefficient::KappaCoefficient(const int two_nu, const int two_L,
                                   const int two_Lp)
    : two_nu(two_nu), two_L(two_L), two_Lp(two_Lp), value(0.) {

  const int nu = two_nu / 2;
  if (nu < 2) {
    throw invalid_argument("nu must be an integer larger than 1.");
  }

  // Avoid division by zero.
  if (!fulfils_triangle_inequality<int>(two_L, two_Lp, two_nu)) {
    value = 0.;
  } else {

    /*
        The expression for \f$\kappa_\nu\f$ given in the literature
        {see, e.g. Eq. (I-7) in \cite FaggHanna1959} use Clebsch-Gordan (CG)
       coefficients, but GSL only provides Wigner-3j symbols. The CG
       coefficients include an additional factor of {see, e.g., Eq. (C.12) in
       \cite Messiah19922}

        \f[
            \left( -1 \right)^{-L + L^\prime - M} \sqrt{2 \nu + 1},
        \f]

        where \f$M\f$ denotes the total magnetic quantum number.
        In the present case, all factors cancel out, and the Wigner-3j symbols
       can be used. Note, however, that \f$M\f$ changes its sign when going from
       the CG coefficient to the Wigner-3j symbol.
    */
    value = -sqrt((double)gsl_sf_fact(nu - 2) / (double)gsl_sf_fact(nu + 2)) *
            gsl_sf_coupling_3j(two_L, two_Lp, two_nu, 2, 2, -4) /
            gsl_sf_coupling_3j(two_L, two_Lp, two_nu, 2, -2, 0);
  }
}

string KappaCoefficient::string_representation(
    const unsigned int n_digits,
    [[maybe_unused]] vector<string> variable_names) const {
  if (n_digits) {
    return float_string_representation(n_digits, value);
  }
  return "\\kappa_{" + to_string(two_nu / 2) + "}" + "\\left(" +
         to_string(two_L / 2) + "," + to_string(two_Lp / 2) + "\\right)";
}
