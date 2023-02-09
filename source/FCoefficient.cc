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

#include <string>

using std::to_string;

#include <gsl/gsl_sf.h>

#include "FCoefficient.hh"
#include "TestUtilities.hh"

FCoefficient::FCoefficient(const int two_nu, const int two_L, const int two_Lp, const int two_j1, const int two_j):
	two_nu(two_nu), two_L(two_L), two_Lp(two_Lp), two_j1(two_j1), two_j(two_j), value(0.){

	double wigner3j{
		gsl_sf_coupling_3j(
			two_L, two_Lp, two_nu,
			2, -2, 0)
	};

	double wigner6j;

	// Shortcut to avoid further calculations.
	if(wigner3j == 0.)
	{
		value = 0.;
	} else {

		wigner6j = gsl_sf_coupling_6j(
				two_j, two_j, two_nu,
				two_Lp, two_L, two_j1);

		// Another shortcut
		if(wigner6j == 0.){
			value =  0.;
		} else {

			value = pow(-1, (two_j1 + two_j)/2 - 1)
			*sqrt((two_L + 1)*(two_Lp + 1)*(two_j + 1)*(two_nu + 1))*wigner3j*wigner6j;

		}
	}
}

bool FCoefficient::is_nonzero(const int two_nu, const int two_L, const int two_Lp, const int two_j1, const int two_j){
	if(
		cg_is_nonzero(
			two_L, two_Lp, two_nu,
			2, -2, 0)
		&& racah_is_nonzero(
			two_j, two_j, two_nu,
			two_Lp, two_L, two_j1
			)
	)
		return true;

	return false;
}

bool FCoefficient::cg_is_nonzero(const int two_j1, const int two_j2, const int two_J, const int two_m1, const int two_m2, const int two_M){

	// Maximum projection of angular momentum.
	if(
		   (abs(two_m1) > two_j1)
		|| (abs(two_m2) > two_j2)
		|| (abs(two_M)  > two_J)
	)
		return false;

	// Conservation of angular momentum for magnetic quantum number.
	if(two_m1 + two_m2 != two_M)
		return false;

	// Triangle inequality for coupling.
	if(!fulfils_triangle_inequality<int>(two_j1, two_j2, two_J))
		return false;

	return true;
}

bool FCoefficient::racah_is_nonzero(const int two_j1, const int two_j2, const int two_j3, const int two_J1, const int two_J2, const int two_J3){

	if(
		   !sum_is_even(two_j1, two_j2, two_j3)
		|| !fulfils_triangle_inequality<int>(two_j1, two_j2, two_j3)
		|| !sum_is_even(two_j1, two_J2, two_J3)
		|| !fulfils_triangle_inequality<int>(two_j1, two_J2, two_J3)
		|| !sum_is_even(two_J1, two_j2, two_J3)
		|| !fulfils_triangle_inequality<int>(two_J1, two_j2, two_J3)
		|| !sum_is_even(two_J1, two_J2, two_j3)
		|| !fulfils_triangle_inequality<int>(two_J1, two_J2, two_j3)
	)
		return false;

	return true;
}

string FCoefficient::string_representation(const unsigned int n_digits, [[maybe_unused]] vector<string> variable_names) const {
	if(n_digits){
		return float_string_representation(n_digits, value);
	}
	string str_rep = "F_{" 
		+ to_string(two_nu/2) 
		+ "}\\left(" 
		+ to_string(two_L/2)
		+ ","
		+ to_string(two_Lp/2)
		+ ",";
	if(two_j1 % 2){
		str_rep += to_string(two_j1)
		+ "/2,"
		+ to_string(two_j)
		+ "/2";
	} else {
		str_rep += to_string(two_j1/2)
		+ ","
		+ to_string(two_j/2);
	}
	str_rep += "\\right)";
	
	return str_rep;
}