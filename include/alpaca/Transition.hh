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
 * \brief Alternate given EMCharacter
 *
 * turns electric into magnetic and vice versa
 *
 * \param em EMCharacter to convert
 *
 * \return Alternate EMCharacter
 */
constexpr EMCharacter alt_character(EMCharacter em) {
  return static_cast<EMCharacter>(-static_cast<short>(em));
}

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
   * \brief Constructor which automatically assigns second transition and does
   * not take information about the EM character
   *
   * The EM characters are initialized as unknown.
   *
   * \param t_L Two times the multipolarity. The second transition corresponds
   * to the next multipolarity order.
   * \param del Multipole mixing ratio.
   * \param del Multipole mixing ratio.
   */
  inline Transition(int t_L, double del = 0.) : Transition(t_L, t_L + 2, del){};

  /**
   * \brief Constructor which automatically assigns second transition
   *
   * \param em Primary EM character. The secondary EM character is assigned
   * automatically.
   * \param t_L Two times the multipolarity. The second
   * transition corresponds to the next multipolarity order.
   * \param del
   * Multipole mixing ratio.
   * \param del Multipole mixing ratio.
   */
  inline Transition(EMCharacter em, int t_L, double del = 0.)
      : Transition(em, t_L, alt_character(em), t_L + 2, del){};

  /**
   * \brief Named constructor for dipole radiation
   *
   * \param del Multipole mixing ratio.
   */
  inline static Transition Dipole(double del = 0.) {
    return Transition(2, del);
  }

  /**
   * \brief Named constructor for E1 radiation
   *
   * \param del Multipole mixing ratio.
   */
  inline static Transition E1(double delta = 0.) {
    return Transition(EMCharacter::electric, 2, delta);
  }

  /**
   * \brief Named constructor for M1 radiation
   *
   * \param del Multipole mixing ratio.
   */
  inline static Transition M1(double delta = 0.) {
    return Transition(EMCharacter::magnetic, 2, delta);
  }

  /**
   * \brief Named constructor for quadrupole radiation
   *
   * \param del Multipole mixing ratio.
   */
  inline static Transition Quadrupole(double delta = 0.) {
    return Transition(4, delta);
  }

  /**
   * \brief Named constructor for E2 radiation
   *
   * \param del Multipole mixing ratio.
   */
  inline static Transition E2(double delta = 0.) {
    return Transition(EMCharacter::electric, 4, delta);
  }

  /**
   * \brief Named constructor for M2 radiation
   *
   * \param del Multipole mixing ratio.
   */
  inline static Transition M2(double delta = 0.) {
    return Transition(EMCharacter::magnetic, 4, delta);
  }

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

inline bool operator==(const Transition &A, const Transition &B) {
  return ((A.em_char == B.em_char) && (A.em_charp == B.em_charp) &&
          (A.two_L == B.two_L) && (A.two_Lp == B.two_Lp) &&
          (A.delta == B.delta));
}

} // namespace alpaca
