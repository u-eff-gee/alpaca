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

#include <sstream>

using std::stringstream;

#include "AlphavCoefficient.hh"

AlphavCoefficient::AlphavCoefficient(const int two_nu, const int two_L, const int two_Lp, const int two_jn, const int two_j):
    two_nu(two_nu), two_L(two_L), two_Lp(two_Lp), two_jn(two_jn), two_j(two_j), constant_f_coefficient(two_nu, two_L, two_L, two_jn, two_j), linear_f_coefficient(two_nu, two_L, two_Lp, two_jn, two_j), quadratic_f_coefficient(two_nu, two_Lp, two_Lp, two_jn, two_j), constant_kappa_coefficient(two_nu, two_L, two_L), linear_kappa_coefficient(two_nu, two_L, two_Lp), quadratic_kappa_coefficient(two_nu, two_Lp, two_Lp)
{
    constant_coefficient = -constant_kappa_coefficient.get_value()*constant_f_coefficient.get_value();
    linear_coefficient = 2.*linear_kappa_coefficient.get_value()*linear_f_coefficient.get_value();
    quadratic_coefficient = quadratic_kappa_coefficient.get_value()*quadratic_f_coefficient.get_value();

}

double AlphavCoefficient::operator()(const double delta) const {
	
	return 
        constant_coefficient + delta*linear_coefficient 
        + delta*delta*quadratic_coefficient;
}

string AlphavCoefficient::string_representation(vector<string> variable_names) const {

    stringstream str_rep;

    str_rep << constant_coefficient 
        << linear_coefficient << " \\times \\delta"
        << quadratic_coefficient << " \\times \\delta^2";

    return str_rep.str();

}