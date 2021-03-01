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

    Copyright (C) 2021 Udo Friman-Gayer
*/

#include <gsl/gsl_sf.h>

#include "EvCoefficient.hh"

EvCoefficient::EvCoefficient(const int two_nu, const EMCharacter em, const int two_L, const EMCharacter emp, const int two_Lp, const int two_jn, const int two_j):
two_nu(two_nu), em(em), two_L(two_L), emp(emp), two_Lp(two_Lp), two_jn(two_jn), two_j(two_j), sign_sigma_L_n((em == magnetic) ? -1 : 1), sign_sigma_Lp_n((emp == magnetic) ? -1 : 1), constant_f_coefficient(two_nu, two_L, two_L, two_jn, two_j), linear_f_coefficient(two_nu, two_L, two_Lp, two_jn, two_j),
quadratic_f_coefficient(two_nu, two_Lp, two_Lp, two_jn, two_j)
{}

double EvCoefficient::operator()(const double delta) const {

    const int nu = two_nu/2;
    const double nu_times_nu_plus_one = nu*(nu+1);
    const int L = two_L/2;
    const double two_L_times_L_plus_one = 2*L*(L+1);
    const int Lp = two_Lp/2;
    const double two_Lp_times_Lp_plus_one = 2*Lp*(Lp+1);

    return (
        sign_sigma_L_n*constant_f_coefficient.get_value()*(nu_times_nu_plus_one*two_L_times_L_plus_one)/(nu_times_nu_plus_one - two_L_times_L_plus_one)
        + 2.*delta*sign_sigma_Lp_n*linear_f_coefficient.get_value()*(Lp - L)*(Lp + L + 1)
        + delta*delta*sign_sigma_Lp_n*quadratic_f_coefficient.get_value()*(nu_times_nu_plus_one*two_Lp_times_Lp_plus_one)/(nu_times_nu_plus_one - two_Lp_times_Lp_plus_one)
    )*gsl_sf_fact(nu-2)/gsl_sf_fact(nu+2);
}
