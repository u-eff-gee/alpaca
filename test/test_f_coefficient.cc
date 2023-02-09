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

#include <cassert>
#include <memory>

using std::make_unique;
using std::unique_ptr;

#include "FCoefficient.hh"
#include "FCoefficientLiteratureValue.hh"
#include "TestUtilities.hh"

vector<FCoefficientLiteratureValue> f_coefficient_values {
	FCoefficientLiteratureValue(2, 0, 2, 2, 4,  0.70710680),
	FCoefficientLiteratureValue(2, 2, 0, 4, 4,  0.00000000),
	FCoefficientLiteratureValue(2, 2, 2, 2, 4, -0.35355340),
	FCoefficientLiteratureValue(2, 2, 2, 4, 4, -1.06066017),
	FCoefficientLiteratureValue(2, 2, 4, 4, 4, -0.35355340),
	FCoefficientLiteratureValue(2, 4, 2, 2, 4,  0.07071068),
	FCoefficientLiteratureValue(2, 4, 2, 4, 4,  0.47434165),
	FCoefficientLiteratureValue(2, 4, 2, 6, 4,  0.84852814),
	FCoefficientLiteratureValue(2, 4, 4, 4, 4,  0.35355340),
	FCoefficientLiteratureValue(2, 4, 4, 6, 4, -0.63245554),
	FCoefficientLiteratureValue(2, 4, 6, 6, 4, -0.42426407),
	FCoefficientLiteratureValue(2, 6, 4, 4, 4, -0.10101526),
	FCoefficientLiteratureValue(2, 6, 4, 6, 4,  0.37796448),
	FCoefficientLiteratureValue(2, 6, 4, 8, 4,  0.95831485),
	FCoefficientLiteratureValue(2, 6, 6, 6, 4,  0.53033009),
	FCoefficientLiteratureValue(2, 6, 6, 8, 4, -0.44821072),
	FCoefficientLiteratureValue(2, 6, 8, 8, 4, -0.42931483),
	FCoefficientLiteratureValue(2, 8, 6, 6, 4, -0.17677670),
	FCoefficientLiteratureValue(2, 8, 6, 8, 4,  0.30618621),
	FCoefficientLiteratureValue(2, 8, 6,10, 4,  0.99999997),
	FCoefficientLiteratureValue(2, 8, 8, 8, 4,  0.60104077),
	FCoefficientLiteratureValue(2, 8, 8,10, 4, -0.34641016),
	FCoefficientLiteratureValue(2, 8,10,10, 4, -0.42426407),
	FCoefficientLiteratureValue(2,10, 8, 8, 4, -0.21856028),
	FCoefficientLiteratureValue(2,10, 8,10, 4,  0.25584086),

	FCoefficientLiteratureValue(5, 7, 12, 12,  4, -0.78963548),
	FCoefficientLiteratureValue(5, 7, 12, 12,  8,  0.29925498),
	FCoefficientLiteratureValue(7, 9,  8,  8, 12, -0.02586818),
	FCoefficientLiteratureValue(7, 9,  8, 10,  4, -0.24618300),
	FCoefficientLiteratureValue(9, 5, 12, 14, 16,  0.45369384),
	FCoefficientLiteratureValue(9, 5, 14, 14,  4, -0.96183098),
	FCoefficientLiteratureValue(17, 1, 16, 16,  4, -1.06646191),
	FCoefficientLiteratureValue(17, 1, 16, 16,  8, 0.95346668),
};

int main(){

	unique_ptr<FCoefficient> f_coef = make_unique<FCoefficient>(4, 0, 0, 0, 0);
	const double epsilon = 1e-7;

	for(auto f : f_coefficient_values){
		f_coef = make_unique<FCoefficient>(f.two_nu, f.two_L, f.two_Lp, f.two_jp, f.two_j);
		// Test whether the numerical value is correct within the given digits of the 
		// literature data.
		test_numerical_equality<double>(f_coef->get_value(), f.value, epsilon);
		// Check whether the prediction that a given F coefficient vanishes is correct.
		assert(FCoefficient::is_nonzero(f.two_nu, f.two_L, f.two_Lp, f.two_jp, f.two_j) == (f_coef->get_value() != 0.));
	}

	f_coef = make_unique<FCoefficient>(4, 2, 2, 0, 2);
	assert(f_coef->string_representation() == "F_{2}\\left(1,1,0,1\\right)");
	assert(f_coef->string_representation(3) == "0.707");

	f_coef = make_unique<FCoefficient>(8, 4, 4, 3, 5);
	assert(f_coef->string_representation() == "F_{4}\\left(2,2,3/2,5/2\\right)");
	assert(f_coef->string_representation(3) == "0.705");

	// Test some of the safety functions that prevent alpaca from calling GSL functions with false
	// sets of arguments.

	// Test Clebsch-Gordan coefficients.
	// This is a correct expression.
	assert(f_coef->cg_is_nonzero(1, 1, 2, 1, -1, 0));
	// Error: m (M) larger than j (J).
	assert(!f_coef->cg_is_nonzero(1, 1, 2, 2, -1, 1));
	assert(!f_coef->cg_is_nonzero(1, 1, 2, 1, -2, -1));
	// The expression below also violates the conservation of angular momentum and the triangle 
	// inequality.
	// However, the magnitudes of J and M are compared first in cg_is_nonzero, so this is truly
	// a test of the condition of interest.
	assert(!f_coef->cg_is_nonzero(1, 1, 2, 1, -1, 3));
	// Error: Conservation of angular momentum violated.
	assert(!f_coef->cg_is_nonzero(1, 1, 2, 1, -1, 2));
	// Error: Triangle inequality violated.
	assert(!f_coef->cg_is_nonzero(1, 1, 3, 1, -1, 0));

	// Test Racah coefficients
	// This is a correct expression.
	assert(f_coef->racah_is_nonzero(1, 1, 2, 1, 1, 2));
	// Error: Sum of j1, j2, j3 not even.
	assert(!f_coef->racah_is_nonzero(1, 1, 1, 1, 1, 2));
	// Error: Sum of j1, J2, J3 not even.
	assert(!f_coef->racah_is_nonzero(1, 1, 1, 1, 1, 2));
	// Error: Sum of J1, j2, J3 not even.
	assert(!f_coef->racah_is_nonzero(1, 1, 2, 2, 1, 1));
	// Error: Sum of J1, J2, J3 not even.
	assert(!f_coef->racah_is_nonzero(1, 1, 2, 1, 1, 1));

	// Error: j1, j2, and j3 do not fulfil the triangle inequality.
	assert(!f_coef->racah_is_nonzero(1, 1, 4, 1, 1, 2));
	// Error: j1, J2, and J3 do not fulfil the triangle inequality.
	assert(!f_coef->racah_is_nonzero(1, 1, 2, 3, 1, 4));
	// Error: J1, j2, and J3 do not fulfil the triangle inequality.
	assert(!f_coef->racah_is_nonzero(3, 1, 2, 1, 1, 4));
	// Error: J1, J2, and J3 do not fulfil the triangle inequality.
	assert(!f_coef->racah_is_nonzero(1, 1, 2, 1, 1, 4));
}
