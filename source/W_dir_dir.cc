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

#include <algorithm>
#include <cmath>
#include <iostream>

#include "State.hh"

#include <gsl/gsl_sf.h>

#include "W_dir_dir.hh"

using std::min;
using std::max;

W_dir_dir::W_dir_dir(const State &ini_sta, const vector<pair<Transition, State>> cas_ste):
W_gamma_gamma(ini_sta, cas_ste)
{
	two_nu_max = calculate_two_nu_max();
	nu_max = two_nu_max/2;
	normalization_factor = calculate_normalization_factor();
	expansion_coefficients = calculate_expansion_coefficients();
}

double W_dir_dir::operator()(const double theta) const {

	double sum_over_nu{0.};

	for(int i = 0; i <= nu_max/2; ++i){
		sum_over_nu += expansion_coefficients[i]*gsl_sf_legendre_Pl(2*i, cos(theta));
	}

	return sum_over_nu*normalization_factor;
}

double W_dir_dir::get_upper_limit() const {

	double upper_limit = 0.;

	for(int i = 0; i <= nu_max/2; ++i){
		upper_limit += fabs(expansion_coefficients[i]);
	}

	return normalization_factor*upper_limit;
}

int W_dir_dir::calculate_two_nu_max() const {
	
	int two_nu_max_Av = calculate_two_nu_max_Av();

	if(n_cascade_steps > 2){
		return min(two_nu_max_Av, calculate_two_nu_max_Uv());
	}
	
	return two_nu_max_Av;
}

int W_dir_dir::calculate_two_nu_max_Av() const {

	return 
	2*min(
		min(
			cascade_steps[0].second.two_J,
			cascade_steps[n_cascade_steps-2].second.two_J
		), 
		min(
			max(cascade_steps[0].first.two_L, cascade_steps[0].first.two_Lp),
			max(cascade_steps[n_cascade_steps-1].first.two_L, cascade_steps[n_cascade_steps-1].first.two_Lp)
		)
	);

}

int W_dir_dir::calculate_two_nu_max_Uv() const {
	
	int min_two_j = cascade_steps[0].second.two_J;

	for(size_t i = 1; i < n_cascade_steps - 1; ++i){
		min_two_j = min(min_two_j, cascade_steps[i].second.two_J);
	}

	return 2*min_two_j;

}

vector<double> W_dir_dir::calculate_expansion_coefficients(){

	vector<double> exp_coef_Av = calculate_expansion_coefficients_Av();

	if(n_cascade_steps > 2){
		vector<double> exp_coef_Uv = calculate_expansion_coefficients_Uv();
		vector<double> exp_coef(exp_coef_Uv.size(), 0.);

		for(size_t i = 0; i < exp_coef_Uv.size(); ++i){
			exp_coef[i] = exp_coef_Av[i]*exp_coef_Uv[i];
		}

		return exp_coef;
	}

	return exp_coef_Av;
}

vector<double> W_dir_dir::calculate_expansion_coefficients_Av(){

	vector<double> exp_coef;
	
	for(int two_nu = 0; two_nu <= two_nu_max; two_nu += 4){
		av_coefficients_excitation.push_back(
			AvCoefficient(two_nu, 
				cascade_steps[0].first.two_L, cascade_steps[0].first.two_Lp,
				initial_state.two_J, cascade_steps[0].second.two_J
			)
		);
		av_coefficients_decay.push_back(
			AvCoefficient(two_nu,
				cascade_steps[n_cascade_steps-1].first.two_L, cascade_steps[n_cascade_steps-1].first.two_Lp,
				cascade_steps[n_cascade_steps-1].second.two_J, cascade_steps[n_cascade_steps-2].second.two_J
			)
		);
		exp_coef.push_back(
			av_coefficients_excitation[two_nu/4](cascade_steps[0].first.delta)
			*av_coefficients_decay[two_nu/4](cascade_steps[n_cascade_steps-1].first.delta)
		);
	}

	return exp_coef;
}

vector<double> W_dir_dir::calculate_expansion_coefficients_Uv(){

	vector<double> exp_coef;
	double uv_coef_product = 1.;

	for(int two_nu = 0; two_nu <= two_nu_max; two_nu += 4){
		uv_coefficients.push_back(vector<UvCoefficient>());
		for(size_t i = 1; i < n_cascade_steps - 1; ++i){
			uv_coefficients[two_nu/4].push_back(
				UvCoefficient(two_nu,
					cascade_steps[i-1].second.two_J,
					cascade_steps[i].first.two_L,
					cascade_steps[i].first.two_Lp,
					cascade_steps[i].first.delta,
					cascade_steps[i].second.two_J
				)
			);
			uv_coef_product = uv_coef_product*uv_coefficients[two_nu/4][i-1].get_value();
		}
		uv_coefficient_products.push_back(uv_coef_product);

		exp_coef.push_back(uv_coef_product);
		uv_coef_product = 1.;
	}

	return exp_coef;
}

double W_dir_dir::calculate_normalization_factor() const {
	
	double norm_fac = 1.;

	for(auto step : cascade_steps){
		norm_fac = norm_fac/(1. + step.first.delta*step.first.delta);
	}

	return norm_fac;

}

string W_dir_dir::string_representation(const unsigned int n_digits, vector<string> variable_names) const {

	const string polar_angle_variable = variable_names.size() ? variable_names[0] : "\\theta";
	const string azimuthal_angle_variable = variable_names.size() ? variable_names[1] : "\\varphi";
	vector<string> delta_variables;
	for(size_t i = 0; i < n_cascade_steps; ++i){
		if(variable_names.size()){
			delta_variables.push_back(variable_names[2+i]);
		} else {
			delta_variables.push_back("\\delta_" + to_string(i+1));
		}
	}

	string str_rep;

	for(int i = 0; i <= nu_max/2; ++i){
		if(i > 0){
			str_rep += "+";
		}
		str_rep += "\\left[" + av_coefficients_excitation[i].string_representation(n_digits, {delta_variables[0]})
			+ "\\right]\\\\";
		if(n_cascade_steps > 2){
			for(size_t j = 0; j < uv_coefficients[i].size(); ++j){
				str_rep += "\\times\\left[" + uv_coefficients[i][j].string_representation(n_digits, {delta_variables[1+j]}) + "\\right]\\\\";
			}
		}
		str_rep += "\\times\\left[" + av_coefficients_decay[i].string_representation(n_digits, {delta_variables[delta_variables.size()-1]})
		+ "\\right]\\\\"
		+ "\\times P_{"
		+ to_string(2*i)
		+ "}\\left[\\cos\\left("
		+ polar_angle_variable
		+ "\\right)\\right]";
		if(i != nu_max/2){
			str_rep += "\\\\";
		}
	}

	return str_rep;
}