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
#include <stdexcept>

#include "alpaca/AngularCorrelation.hh"
#include "alpaca/State.hh"
#include "alpaca/Transition.hh"

using alpaca::AngularCorrelation;
using alpaca::EMCharacter;
using alpaca::Parity;
using alpaca::State;
using alpaca::Transition;

/**
 * \brief Test the input validation of the AngularCorrelation class.
 *
 * Test by calling the constructor with various invalid inputs.
 * The basic example is a
 *
 * 0^+ -> 1^+ -> 0^+
 *
 * cascade.
 */
int main() {

  bool error_thrown = false;

  // Error: First EM character not given
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: First EM character not given for second transition
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Second EM character not given
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Second EM character not given for second transition
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: EM character given, but first parity missing
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: EM character given, but second parity missing
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: EM character given, but parity missing for final state of second
  // transition
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::unknown)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Triangle inequality violated for first transition
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Triangle inequality violated for second transition
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Not an error: Triangle inequality fulfilled by second transition
  {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 10, EMCharacter::electric, 2, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  }

  // Error: First electromagnetic character wrong
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::magnetic, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: First electromagnetic character wrong for second transition
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::magnetic, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Second electromagnetic character wrong
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Second electromagnetic character wrong for second transition
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(2, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::electric, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Mixing of half-integer and integer spins
  try {
    AngularCorrelation ang_corr(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(1, Parity::negative)},
         {Transition(EMCharacter::electric, 2, EMCharacter::magnetic, 4, 0.),
          State(0, Parity::positive)}});
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Simplified constructor. Check transition inference with various
  // possibilities. Error: Transition between spin-0 states
  try {
    AngularCorrelation ang_corr_0_1_0{
        State(0, Parity::unknown),
        {State(0, Parity::negative), State(0, Parity::unknown)}};
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Too few steps in cascade
  try {
    AngularCorrelation ang_corr_0_1{State(0, Parity::unknown),
                                    {
                                        State(2, Parity::negative),
                                    }};
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  // Error: Same as above with more general constructor.
  try {
    AngularCorrelation ang_corr_0p_1p{
        State(0, Parity::positive),
        {
            {Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
             State(2, Parity::positive)},
        }};
  } catch (const std::invalid_argument &e) {
    error_thrown = true;
  }

  assert(error_thrown);
  error_thrown = false;

  AngularCorrelation ang_corr_0_1_0{
      State(0, Parity::positive),
      {State(2, Parity::negative), State(6, Parity::negative)}};

  vector<pair<Transition, State>> cas_ste = ang_corr_0_1_0.get_cascade_steps();

  assert(cas_ste[0].first.two_L == 2);
  assert(cas_ste[0].first.em_char == EMCharacter::electric);
  assert(cas_ste[0].first.two_Lp == 4);
  assert(cas_ste[0].first.em_charp == EMCharacter::magnetic);
  assert(cas_ste[0].first.delta == 0.);
  assert(cas_ste[1].first.two_L == 4);
  assert(cas_ste[1].first.em_char == EMCharacter::electric);
  assert(cas_ste[1].first.two_Lp == 6);
  assert(cas_ste[1].first.em_charp == EMCharacter::magnetic);
  assert(cas_ste[1].first.delta == 0.);

  AngularCorrelation ang_corr_3_5_3_9 =
      AngularCorrelation{State(3, Parity::positive),
                         {State(5, Parity::unknown), State(3, Parity::negative),
                          State(9, Parity::negative)}};

  cas_ste = ang_corr_3_5_3_9.get_cascade_steps();

  assert(cas_ste[0].first.em_char == EMCharacter::unknown);
  assert(cas_ste[0].first.em_charp == EMCharacter::unknown);
  assert(cas_ste[1].first.em_char == EMCharacter::unknown);
  assert(cas_ste[1].first.em_charp == EMCharacter::unknown);
  assert(cas_ste[2].first.em_char == EMCharacter::magnetic);
  assert(cas_ste[2].first.em_charp == EMCharacter::electric);

  // Check that the correct transition is inferred between states of equal spin.
  // The triangle inequality gives a monopole transition as the lowest possible
  // multipolarity in that case, but this is not possible for a single photon.
  AngularCorrelation ang_corr_0_4_4 =
      AngularCorrelation{State(0, Parity::positive),
                         {
                             State(4, Parity::positive),
                             State(4, Parity::positive),
                         }};

  cas_ste = ang_corr_0_4_4.get_cascade_steps();

  assert(cas_ste[0].first.em_char == EMCharacter::electric);
  assert(cas_ste[0].first.two_L == 4);
  assert(cas_ste[0].first.em_charp == EMCharacter::magnetic);
  assert(cas_ste[0].first.two_Lp == 6);
  assert(cas_ste[1].first.em_char == EMCharacter::magnetic);
  assert(cas_ste[1].first.two_L == 2);
  assert(cas_ste[1].first.em_charp == EMCharacter::electric);
  assert(cas_ste[1].first.two_Lp == 4);

  AngularCorrelation ang_corr_0_4m_4 =
      AngularCorrelation{State(0, Parity::positive),
                         {
                             State(4, Parity::negative),
                             State(4, Parity::positive),
                         }};

  cas_ste = ang_corr_0_4m_4.get_cascade_steps();

  assert(cas_ste[0].first.em_char == EMCharacter::magnetic);
  assert(cas_ste[0].first.two_L == 4);
  assert(cas_ste[0].first.em_charp == EMCharacter::electric);
  assert(cas_ste[0].first.two_Lp == 6);
  assert(cas_ste[1].first.em_char == EMCharacter::electric);
  assert(cas_ste[1].first.two_L == 2);
  assert(cas_ste[1].first.em_charp == EMCharacter::magnetic);
  assert(cas_ste[1].first.two_Lp == 4);
}
