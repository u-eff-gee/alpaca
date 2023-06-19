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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include "alpaca/W_pol_dir.hh"

using std::string;
using std::to_string;

namespace alpaca {

W_pol_dir::W_pol_dir(const State &ini_sta,
                     const vector<pair<Transition, State>> &cas_ste)
    : W_gamma_gamma(ini_sta, cas_ste), w_dir_dir(W_dir_dir(ini_sta, cas_ste)) {

  two_nu_max = w_dir_dir.get_two_nu_max();
  nu_max = two_nu_max / 2;
  expansion_coefficients = calculate_expansion_coefficients();
  normalization_factor = w_dir_dir.get_normalization_factor();
}

double W_pol_dir::operator()(const double theta, const double phi) const {

  double sum_over_nu{0.};

  for (int i = 1; i <= nu_max / 2; ++i) {
    sum_over_nu += expansion_coefficients[static_cast<size_t>(i - 1)] *
                   gsl_sf_legendre_Plm(2 * i, 2, cos(theta));
  }

  int polarization_sign = 1;
  if (cascade_steps[0].first.em_charp == EMCharacter::magnetic) {
    polarization_sign = -1;
  }

  return w_dir_dir(theta) + polarization_sign * cos(2. * phi) * sum_over_nu *
                                w_dir_dir.get_normalization_factor();
}

double W_pol_dir::get_upper_limit() const {

  double upper_limit = 0.;

  const double associated_Legendre_upper_limit_factor = 4. * pow(M_1_PI, 0.75);

  for (int i = 1; i <= nu_max / 2; ++i) {
    upper_limit += fabs(expansion_coefficients[static_cast<size_t>(i - 1)]) *
                   associated_Legendre_upper_limit_factor *
                   sqrt(gsl_sf_fact(static_cast<unsigned int>(2 * i + 2)) /
                        gsl_sf_fact(static_cast<unsigned int>(2 * i - 2)));
  }

  return w_dir_dir.get_upper_limit() +
         upper_limit * w_dir_dir.get_normalization_factor();
}

vector<double> W_pol_dir::calculate_expansion_coefficients() {

  vector<double> exp_coef_alphav_Av =
      calculate_expansion_coefficients_alphav_Av();

  if (n_cascade_steps > 2) {
    vector<double> exp_coef_Uv = w_dir_dir.get_Uv_coefficient_products();
    vector<double> exp_coef(exp_coef_Uv.size(), 0.);

    for (size_t i = 1; i < exp_coef_Uv.size(); ++i) {
      exp_coef[i - 1] = exp_coef_alphav_Av[i - 1] * exp_coef_Uv[i];
    }

    return exp_coef;
  }

  return exp_coef_alphav_Av;
}

vector<double> W_pol_dir::calculate_expansion_coefficients_alphav_Av() {

  vector<double> exp_coef;

  for (int two_nu = 4; two_nu <= two_nu_max; two_nu += 4) {
    alphav_coefficients.emplace_back(
        two_nu, cascade_steps[0].first.two_L, cascade_steps[0].first.two_Lp,
        initial_state.two_J, cascade_steps[0].second.two_J);
    av_coefficients.emplace_back(
        two_nu, cascade_steps[n_cascade_steps - 1].first.two_L,
        cascade_steps[n_cascade_steps - 1].first.two_Lp,
        cascade_steps[n_cascade_steps - 1].second.two_J,
        cascade_steps[n_cascade_steps - 2].second.two_J);
    exp_coef.push_back(alphav_coefficients.at(static_cast<size_t>(
                           two_nu / 4 - 1))(cascade_steps[0].first.delta) *
                       av_coefficients.at(static_cast<size_t>(two_nu / 4 - 1))(
                           cascade_steps[n_cascade_steps - 1].first.delta));
  }

  return exp_coef;
}

string W_pol_dir::string_representation(const int n_digits,
                                        vector<string> variable_names) const {

  const string polar_angle_variable =
      variable_names.empty() ? "\\theta" : variable_names[0];
  const string azimuthal_angle_variable =
      variable_names.empty() ? "\\varphi" : variable_names[1];
  vector<string> delta_variables;
  for (size_t i = 0; i < n_cascade_steps; ++i) {
    if (variable_names.empty()) {
      delta_variables.push_back("\\delta_" + to_string(i + 1));
    } else {
      delta_variables.push_back(variable_names[2 + i]);
    }
  }

  const vector<vector<UvCoefficient>> uv_coefficients =
      w_dir_dir.get_Uv_coefficients();

  string str_rep =
      w_dir_dir.string_representation(n_digits, variable_names) + "\\\\";
  str_rep +=
      cascade_steps[0].first.em_charp == EMCharacter::magnetic ? "+" : "-";
  str_rep += R"(\cos\left(2)" + azimuthal_angle_variable +
             R"(\right)\left\{\right.\\)";

  for (int i = 1; i <= nu_max / 2; ++i) {
    if (i > 1) {
      str_rep += "+";
    }

    str_rep +=
        R"(\left[)" +
        alphav_coefficients[static_cast<size_t>(i - 1)].string_representation(
            n_digits, {delta_variables[0]}) +
        R"(\right]\\)";
    if (n_cascade_steps > 2) {
      for (size_t j = 0; j < uv_coefficients[static_cast<size_t>(i)].size();
           ++j) {
        str_rep +=
            R"(\times\left[)" +
            uv_coefficients[static_cast<size_t>(i)][j].string_representation(
                n_digits, {delta_variables[1 + j]}) +
            R"(\right]\\)";
      }
    }
    str_rep +=
        R"(\times\left[)" +
        av_coefficients[static_cast<size_t>(i - 1)].string_representation(
            n_digits, {delta_variables[delta_variables.size() - 1]}) +
        R"(\right]\\\times P_{)" + to_string(2 * i) +
        R"(}^{\left|2\right|}\left[\cos\left()" + polar_angle_variable +
        "\\right)\\right]";
    if (i != nu_max / 2) {
      str_rep += "\\\\";
    }
  }
  str_rep += R"(\left.\right\})";

  return str_rep;
}

} // namespace alpaca
