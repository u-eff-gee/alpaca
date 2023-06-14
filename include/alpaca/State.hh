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

using std::invalid_argument;
using std::string;

namespace alpaca {

/**
 * \brief Enum for the possible values of the parity quantum number.
 */
enum class Parity : short { negative = -1, positive = 1, unknown = 0 };

/**
 * \brief Struct to store properties of a nuclear state.
 *
 * The state is characterized by its angular momentum, parity quantum number,
 * and excitation energy. The excitation energy is assumed to be defined with
 * respect to the ground state. This means that the ground state has an energy
 * of 0 MeV, and all excited states must have a positive nonzero energy.
 */
struct State {
  /**
   * \brief Constructor which does not require energy or parity information
   *
   * The parity quantum number is initialized as unknown.
   * The excitation energy is initialized as 0 MeV.
   *
   * \param t_J Two times the angular momentum quantum number in units of the
   * reduced Planck constant.
   *
   * \throw std::invalid_argument if two_J is invalid
   */
  explicit State(const int t_J)
      : two_J(check_two_J(t_J)), parity(Parity::unknown),
        excitation_energy(0.){};
  State(double t_J) = delete;

  /**
   * \brief Constructor which does not require energy information
   *
   * The excitation energy is initialized as 0 MeV.
   *
   * \param t_J Two times the angular momentum quantum number in units of the
   * reduced Planck constant. \param parity Parity quantum number.
   *
   * \throw std::invalid_argument if two_J is invalid
   */
  State(const int t_J, const Parity p)
      : two_J(check_two_J(t_J)), parity(p), excitation_energy(0.){};
  State(double t_J, Parity p) = delete;
  /**
   * \brief Constructor which does not require parity information
   *
   * The parity quantum number is initialized as unknown.
   *
   * \param t_J Two times the angular momentum quantum number in units of the
   * reduced Planck constant. \param excitation_energy Excitation energy of the
   * state with respect to the ground state in MeV.
   *
   * \throw std::invalid_argument if two_J is invalid
   */
  State(const int t_J, const Parity p, const double e_x)
      : two_J(check_two_J(t_J)), parity(p),
        excitation_energy(check_excitation_energy(e_x)){};
  State(double t_J, Parity p, double e_x) = delete;
  /**
   * \brief Constructor
   *
   * \param t_J Two times the angular momentum quantum number in units of the
   * reduced Planck constant. \param parity Parity quantum number. \param
   * excitation_energy Excitation energy of the state with respect to the ground
   * state in MeV.
   *
   * \throw std::invalid_argument if two_J is invalid
   */
  State(const int t_J, const double e_x)
      : two_J(check_two_J(t_J)), parity(Parity::unknown),
        excitation_energy(check_excitation_energy(e_x)){};
  State(double t_J, double e_x) = delete;

  int two_J; /**< Two times the angular momentum quantum number in units of the
                reduced Planck constant. */
  Parity parity;            /**< Parity quantum number. */
  double excitation_energy; /**< Excitation energy of the state with respect to
                               the ground state in MeV. */

  /**
   * \brief String representation of parities.
   *
   * \return "+" or "-"
   *
   * \throw runtime_error if parity is neither positive nor negative.
   */
  [[nodiscard]] inline string parity_str_rep() const {
    switch (parity) {
    case Parity::negative:
      return "-";
    case Parity::positive:
      return "+";
    default:
      throw std::runtime_error("No string representation for unknown parity.");
    }
  }

  /**
   * \brief String representation of angular momentum quantum numbers.
   *
   * \return String representation
   */
  [[nodiscard]] inline string spin_str_rep() const {
    if (two_J % 2 == 0) {
      return std::to_string(two_J / 2);
    }

    return std::to_string(two_J) + "/2";
  }

  /**
   * \brief String representation of a State.
   *
   * If the parity is unknown, it is omitted.
   *
   * \param state State
   *
   * \return String representation
   */
  [[nodiscard]] string str_rep() const {
    if (parity != Parity::unknown) {
      return spin_str_rep() + "^" + parity_str_rep();
    }

    return spin_str_rep();
  }

  /**
   * \brief Ensure that given angular momentum quantum number is valid.
   *
   * The reason why two_J was defined as an 'int' and not an 'unsigned int' is
   * because the GSL \cite Galassi2009 functions accept 'int'.
   *
   * \param int two_J
   *
   * \returns two_J, if it is valid
   *
   * \throw std::invalid_argument if two_J is invalid
   */
  constexpr static int check_two_J(const int two_J) {
    if (two_J < 0) {
      throw invalid_argument("two_J must be a nonnegative integer.");
    }

    return two_J;
  }

  /**
   * \brief Ensure that given excitation energy is valid.
   *
   * \param double e_x
   *
   * \returns e_x, if it is valid
   *
   * \throw std::invalid_argument if e_x is invalid.
   */
  constexpr static double check_excitation_energy(const double e_x) {
    if (e_x < 0.) {
      throw std::invalid_argument("Excitation energy must not be negative.");
    }

    return e_x;
  }

  friend bool operator==(State const &, State const &) = default;
};

} // namespace alpaca
