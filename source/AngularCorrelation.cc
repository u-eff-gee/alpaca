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

#include <stdexcept>
#include <string>

using std::invalid_argument;
using std::to_string;

#include "alpaca/AngularCorrelation.hh"
#include "alpaca/EulerAngleRotation.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/W_dir_dir.hh"
#include "alpaca/W_pol_dir.hh"

namespace alpaca {

AngularCorrelation::AngularCorrelation(
    const State ini_sta, const vector<pair<Transition, State>> cas_ste)
    : w_gamma_gamma(nullptr) {
  check_cascade(ini_sta, cas_ste);

  if (cas_ste[0].first.em_char == EMCharacter::unknown) {
    w_gamma_gamma = std::make_unique<W_dir_dir>(ini_sta, cas_ste);
  } else {
    w_gamma_gamma = std::make_unique<W_pol_dir>(ini_sta, cas_ste);
  }
}

AngularCorrelation::AngularCorrelation(const State ini_sta,
                                       const vector<State> cas_sta)
    : w_gamma_gamma(nullptr) {

  if (cas_sta.size() < 2) {
    throw invalid_argument(
        "Cascade must have at least two transition - state pairs.");
  }

  vector<pair<Transition, State>> cascade_steps;

  cascade_steps.emplace_back(infer_transition({ini_sta, cas_sta[0]}),
                             cas_sta[0]);

  for (size_t i = 0; i < cas_sta.size() - 1; ++i) {
    cascade_steps.emplace_back(infer_transition({cas_sta[i], cas_sta[i + 1]}),
                               cas_sta[i + 1]);
  }

  if (cascade_steps[0].first.em_char == EMCharacter::unknown) {
    w_gamma_gamma = std::make_unique<W_dir_dir>(ini_sta, cascade_steps);
  } else {
    w_gamma_gamma = std::make_unique<W_pol_dir>(ini_sta, cascade_steps);
  }
}

Transition
AngularCorrelation::infer_transition(const pair<State, State> states) const {

  if (states.first.two_J == 0 and states.second.two_J == 0) {
    throw invalid_argument(
        "An electromagnetic transition between two spin-0 states with the "
        "absorption/emission of a single photon is not possible.");
  }

  int two_L = abs(states.first.two_J - states.second.two_J);
  if (two_L == 0) {
    two_L += 2;
  }
  EMCharacter em = EMCharacter::unknown;
  EMCharacter emp = EMCharacter::unknown;
  if (states.first.parity != Parity::unknown &&
      states.second.parity != Parity::unknown) {
    if (two_L % 4 == 0) {
      if (states.first.parity == states.second.parity) {
        em = EMCharacter::electric;
        emp = EMCharacter::magnetic;
      } else {
        em = EMCharacter::magnetic;
        emp = EMCharacter::electric;
      }
    } else {
      if (states.first.parity == states.second.parity) {
        em = EMCharacter::magnetic;
        emp = EMCharacter::electric;
      } else {
        em = EMCharacter::electric;
        emp = EMCharacter::magnetic;
      }
    }
  }

  return Transition(em, two_L, emp, two_L + 2, 0.);
}

double AngularCorrelation::operator()(const double theta,
                                      const double phi) const {
  return w_gamma_gamma->operator()(theta, phi);
}

void AngularCorrelation::check_cascade(
    const State ini_sta, const vector<pair<Transition, State>> cas_ste) const {

  if (cas_ste.size() < 2) {
    throw invalid_argument(
        "Cascade must have at least two transition - state pairs.");
  }

  check_angular_momenta(ini_sta, cas_ste);

  check_triangle_inequalities(ini_sta, cas_ste);

  check_em_transitions(ini_sta, cas_ste);
}

void AngularCorrelation::check_angular_momenta(
    const State ini_sta, const vector<pair<Transition, State>> cas_ste) const {

  const int even_odd = ini_sta.two_J % 2;

  for (const auto &step : cas_ste) {
    if (step.second.two_J % 2 != even_odd) {
      throw invalid_argument(
          "Unphysical mixing of half-integer and integer spins in cascade.");
    }
  }
}

void AngularCorrelation::check_triangle_inequalities(
    const State ini_sta, const vector<pair<Transition, State>> cas_ste) const {

  if (!fulfils_triangle_inequality<int>(ini_sta.two_J, cas_ste[0].second.two_J,
                                        cas_ste[0].first.two_L) &&
      !fulfils_triangle_inequality<int>(ini_sta.two_J, cas_ste[0].second.two_J,
                                        cas_ste[0].first.two_Lp)) {
    throw invalid_argument(
        "Triangle inequality selection rule not fulfilled for any "
        "multipolarity of transition #1: " +
        cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
  }

  for (size_t i = 1; i < cas_ste.size(); ++i) {
    if (!fulfils_triangle_inequality<int>(cas_ste[i - 1].second.two_J,
                                          cas_ste[i].second.two_J,
                                          cas_ste[i].first.two_L) &&
        !fulfils_triangle_inequality<int>(cas_ste[i - 1].second.two_J,
                                          cas_ste[i].second.two_J,
                                          cas_ste[i].first.two_Lp)) {
      throw invalid_argument(
          "Triangle inequality selection rule not fulfilled for any "
          "multipolarity of transition #" +
          cas_ste[i].first.str_rep(cas_ste[i - 1].second, cas_ste[i].second));
    }
  }
}

void AngularCorrelation::check_em_transitions(
    const State ini_sta, const vector<pair<Transition, State>> cas_ste) const {

  if (ini_sta.parity != Parity::unknown &&
      cas_ste[0].second.parity != Parity::unknown) {
    if (cas_ste[0].first.em_char != EMCharacter::unknown) {
      if (!valid_em_character(ini_sta.parity, cas_ste[0].second.parity,
                              cas_ste[0].first.two_L,
                              cas_ste[0].first.em_char)) {
        throw invalid_argument(
            "Incorrect electromagnetic character '" +
            Transition::em_str_rep(cas_ste[0].first.em_char) +
            "' for transition #1: " +
            cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
      }

      if (cas_ste[0].first.em_charp != EMCharacter::unknown) {
        if (!valid_em_character(ini_sta.parity, cas_ste[0].second.parity,
                                cas_ste[0].first.two_Lp,
                                cas_ste[0].first.em_charp)) {
          throw invalid_argument(
              "Incorrect electromagnetic character '" +
              Transition::em_str_rep(cas_ste[0].first.em_charp) +
              "' for transition #1: " +
              cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
        }
      } else {
        throw invalid_argument(
            "Only one electromagnetic character defined for transition #1: " +
            cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
      }
    }

    if (cas_ste[0].first.em_char == EMCharacter::unknown &&
        cas_ste[0].first.em_charp != EMCharacter::unknown) {
      throw invalid_argument(
          "Only one electromagnetic character defined for transition #1: " +
          cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
    }
  } else if (cas_ste[0].first.em_char != EMCharacter::unknown ||
             cas_ste[0].first.em_charp != EMCharacter::unknown) {
    throw invalid_argument(
        "Electromagnetic character defined, but one or both parities missing "
        "for transition #1: " +
        cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
  }

  for (size_t i = 1; i < cas_ste.size(); ++i) {
    if (cas_ste[i - 1].second.parity != Parity::unknown &&
        cas_ste[i].second.parity != Parity::unknown) {
      if (cas_ste[i].first.em_char != EMCharacter::unknown) {
        if (!valid_em_character(
                cas_ste[i - 1].second.parity, cas_ste[i].second.parity,
                cas_ste[i].first.two_L, cas_ste[i].first.em_char)) {
          throw invalid_argument(
              "Incorrect electromagnetic character '" +
              Transition::em_str_rep(cas_ste[0].first.em_char) +
              "' for transition #" + to_string(i + 1) + ": " +
              cas_ste[i].first.str_rep(cas_ste[i - 1].second,
                                       cas_ste[i].second));
        }
      }

      if (cas_ste[i].first.em_charp != EMCharacter::unknown) {
        if (!valid_em_character(
                cas_ste[i - 1].second.parity, cas_ste[i].second.parity,
                cas_ste[i].first.two_Lp, cas_ste[i].first.em_charp)) {
          throw invalid_argument(
              "Incorrect electromagnetic character '" +
              Transition::em_str_rep(cas_ste[0].first.em_charp) +
              "' for transition #" + to_string(i + 1) + ": " +
              cas_ste[i].first.str_rep(cas_ste[i - 1].second,
                                       cas_ste[i].second));
        }
      } else {
        throw invalid_argument(
            "Only one electromagnetic character defined for transition #" +
            to_string(i + 1) + ": " +
            cas_ste[i].first.str_rep(cas_ste[i - 1].second, cas_ste[i].second));
      }

      if (cas_ste[i].first.em_char == EMCharacter::unknown &&
          cas_ste[i].first.em_charp != EMCharacter::unknown) {
        throw invalid_argument(
            "Only one electromagnetic character defined for transition #" +
            to_string(i + 1) + ": " +
            cas_ste[i].first.str_rep(cas_ste[i - 1].second, cas_ste[i].second));
      }
    } else if (cas_ste[i].first.em_char != EMCharacter::unknown ||
               cas_ste[i].first.em_charp != EMCharacter::unknown) {
      throw invalid_argument(
          "Electromagnetic character defined, but one or both parities missing "
          "for transition #" +
          to_string(i + 1) + ": " +
          cas_ste[i].first.str_rep(cas_ste[i - 1].second, cas_ste[i].second));
    }
  }
}

bool AngularCorrelation::valid_em_character(const Parity p0, const Parity p1,
                                            const int two_L,
                                            const EMCharacter em) const {

  if (p0 == p1) {
    if ((two_L / 2) % 2 == 0) {
      if (em != EMCharacter::electric) {
        return false;
      }
    } else if (em != EMCharacter::magnetic) {
      return false;
    }

    return true;
  }

  if ((two_L / 2) % 2 == 0) {
    if (em != EMCharacter::magnetic) {
      return false;
    }
  } else if (em != EMCharacter::electric) {
    return false;
  }

  return true;
}

} // namespace alpaca
