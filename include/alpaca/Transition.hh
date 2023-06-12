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

#include <stdexcept>
#include <string>

#include "alpaca/State.hh"

using std::invalid_argument;
using std::runtime_error;
using std::string;

namespace alpaca {

/**
 * \brief Enum for the possible values of the electromagnetic (EM) character.
 */
enum class EMCharacter : short { electric = -1, magnetic = 1, unknown = 0 };

/**
 * \brief Struct to store properties of an EM transition between nuclear states.
 *
 * The transition can have two different multipolarities with their associated
 * EM character, whose relative intensity is given by the multipole mixing
 * ratio.
 */
struct Transition {
  /**
   * \brief Constructor which does not take information about the EM character
   *
   * The EM characters are initialized as unknown.
   *
   * \param t_L Two times the multipolarity.
   * \param t_Lp Two times the alternative multipolarity.
   * Must be different from t_L.
   * \param del Multipole mixing ratio.
   *
   * \throw invalid_argument if an invalid value for \f$2 L\f$ or \f$2
   * L^\prime\f$ was given, or if the two are equal.
   */
  Transition(const int t_L, const int t_Lp, const double del = 0.);
  /**
   * \brief Constructor
   *
   * \param em Primary EM character.
   * \param t_L Two times the primary multipolarity.
   * \param emp Secondary EM character.
   * \param t_Lp Two times the secondary multipolarity.
   * Must be different from t_L.
   * \param del Multipole mixing ratio.
   *
   * \throw invalid_argument if an invalid value for \f$2 L\f$ or \f$2
   * L^\prime\f$ was given, or if the two are equal.
   */
  Transition(const EMCharacter em, const int t_L, const EMCharacter emp,
             const int t_Lp, const double del);

  /**
   * \brief String representation of EM characters.
   *
   * \param em \f$\lambda\f$, EM character
   *
   * \return "E" or "M"
   *
   * \throw runtime_error if em is neither electric nor magnetic.
   */
  static string em_str_rep(const EMCharacter em) {

    if (em == EMCharacter::electric) {
      return "E";
    }
    if (em == EMCharacter::magnetic) {
      return "M";
    }

    throw runtime_error(
        "No string representation for unknown electromagnetic character.");
  }

  /**
   * \brief String representation of a transition between two states.
   *
   * If parities or EM characters are unknown, they are omitted.
   * At the moment, the secondary multipolarity will be shown even if the
   * transition is pure.
   *
   * \param initial_state Initial state of the transition
   * \param final_state Final state of the transition
   *
   * \return String representation
   */
  string str_rep(const State initial_state, const State final_state) const;

  EMCharacter em_char;  /**< Primary EM character. */
  int two_L;            /**< Two times the primary multipolarity. */
  EMCharacter em_charp; /**< Secondary EM character. */
  int two_Lp;           /**< Two times the secondary multipolarity. */
  double delta;         /**< Multipole mixing ratio. */

  /**
   * \brief Ensure that given multipolarity is valid.
   *
   * The reason why two_L was defined as an 'int' and not an 'unsigned int' is
   * because the GSL \cite Galassi2009 functions accept 'int'.
   *
   * \param int two_L
   *
   * \returns two_L, if it is valid
   *
   * \throw std::invalid_argument if two_L is invalid
   */
  int check_two_L(const int two_L) const;
};

} // namespace alpaca
