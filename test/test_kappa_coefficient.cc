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

#include <cassert>
#include <memory>

using std::make_unique;
using std::unique_ptr;

#include <stdexcept>

#include "KappaCoefficient.hh"
#include "KappaCoefficientLiteratureValue.hh"
#include "TestUtilities.hh"

/*
	From Tab. II(b) of Ref. \cite FaggHanna1959
*/
vector<KappaCoefficientLiteratureValue> kappa_coefficient_values {
	KappaCoefficientLiteratureValue(4, 2, 2, -1./2.),
	KappaCoefficientLiteratureValue(4, 2, 4, -1./6.),
	KappaCoefficientLiteratureValue(4, 2, 6, -1./12.),
	KappaCoefficientLiteratureValue(4, 4, 4, 1./2.),
	KappaCoefficientLiteratureValue(4, 4, 6, -1./4.),
	KappaCoefficientLiteratureValue(4, 6, 6, 1./3.),

	KappaCoefficientLiteratureValue(6, 2, 4, -1./6.),
	KappaCoefficientLiteratureValue(6, 2, 6, -1./12.),
	KappaCoefficientLiteratureValue(6, 4, 4, 0.),
	KappaCoefficientLiteratureValue(6, 4, 6, 1./4.),
	KappaCoefficientLiteratureValue(6, 6, 6, 0.),

	KappaCoefficientLiteratureValue(8, 2, 6, -1./12.),
	KappaCoefficientLiteratureValue(8, 4, 4, -1./12.),
	KappaCoefficientLiteratureValue(8, 4, 6, -1./60.),
	KappaCoefficientLiteratureValue(8, 6, 6, 1./3.),

	KappaCoefficientLiteratureValue(10, 4, 6, -1./20.),
	KappaCoefficientLiteratureValue(10, 6, 6, 0.),

	KappaCoefficientLiteratureValue(12, 6, 6, -1./30.),
};

int main(){

	unique_ptr<KappaCoefficient> kappa_coef = make_unique<KappaCoefficient>(4, 0, 0);
	const double epsilon = 1e-7;

	for(auto k : kappa_coefficient_values){
		kappa_coef = make_unique<KappaCoefficient>(k.two_nu, k.two_L, k.two_Lp);
		// Test whether the numerical value is correct within the given digits of the 
		// literature data.
		test_numerical_equality<double>(kappa_coef->get_value(), k.value, epsilon);
	}

	[[maybe_unused]] bool error_thrown{false};
	try{
		kappa_coef = make_unique<KappaCoefficient>(0, 2, 2);
	} catch(const std::invalid_argument &e){
		error_thrown = true;
	}
	assert(error_thrown);

	kappa_coef = make_unique<KappaCoefficient>(4, 2, 2);
	assert(kappa_coef->string_representation() == "\\kappa_{2}\\left(1,1\\right)");
	assert(kappa_coef->string_representation(3) == "\\left(-0.5\\right)");
}
