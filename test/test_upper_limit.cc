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

#include <array>
#include <cassert>
#include <vector>

using std::array;
using std::vector;

#include <gsl/gsl_sf.h>

#include "alpaca/SpherePointSampler.hh"
#include "alpaca/State.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/Transition.hh"
#include "alpaca/W_dir_dir.hh"
#include "alpaca/W_gamma_gamma.hh"
#include "alpaca/W_pol_dir.hh"

using alpaca::EMCharacter;
using alpaca::Parity;
using alpaca::SpherePointSampler;
using alpaca::State;
using alpaca::Transition;
using alpaca::W_dir_dir;
using alpaca::W_gamma_gamma;
using alpaca::W_pol_dir;

vector<W_gamma_gamma *> ang_corrs{
    // // 0 -> 1 -> 0
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::unknown)}}),

    // 0 -> 2 -> 0
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(0, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(0, Parity::unknown)}}),

    // 0 -> 1 -> 1
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(2, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(2, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(2, Parity::unknown)}}),

    // 0 -> 1 -> 2
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(4, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(4, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(4, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(4, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(4, Parity::unknown)}}),

    // 0 -> 1 -> 3
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(6, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(6, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 1.),
          State(6, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, -1.),
          State(6, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 100.),
          State(6, Parity::unknown)}}),

    // 0 -> 2 -> 1
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(2, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(2, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(2, Parity::unknown)}}),

    // 0 -> 2 -> 2
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(4, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(4, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(4, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(4, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(4, Parity::unknown)}}),

    // 0 -> 2 -> 3
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(6, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(6, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(6, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(6, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(6, Parity::unknown)}}),

    // 0 -> 2 -> 4
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(8, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(8, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 1.),
          State(8, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, -1.),
          State(8, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 4, EMCharacter::electric, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 100.),
          State(8, Parity::unknown)}}),

    // 1/2 -> 5/2 -> 1/2
    new W_dir_dir(
        State(1, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(5, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(1, Parity::unknown)}}),
    new W_pol_dir(
        State(1, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(1, Parity::unknown)}}),
    new W_pol_dir(
        State(1, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 1.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 1.),
          State(1, Parity::unknown)}}),
    new W_pol_dir(
        State(1, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, -1.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, -1.),
          State(1, Parity::unknown)}}),
    new W_pol_dir(
        State(1, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 100.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 100.),
          State(1, Parity::unknown)}}),

    // 3/2 -> 3/2 -> 3/2
    new W_dir_dir(
        State(3, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(3, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(3, Parity::unknown)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(3, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(3, Parity::positive)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 1.),
          State(3, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(3, Parity::positive)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, -1.),
          State(3, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(3, Parity::positive)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 100.),
          State(3, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(3, Parity::positive)}}),

    // 3/2 -> 5/2 -> 3/2
    new W_dir_dir(
        State(3, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(5, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(3, Parity::unknown)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(3, Parity::positive)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 1.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(3, Parity::positive)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, -1.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(3, Parity::positive)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 100.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(3, Parity::positive)}}),

    // 3/2 -> 7/2 -> 3/2
    new W_dir_dir(
        State(3, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(7, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(3, Parity::unknown)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(7, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(3, Parity::unknown)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 1.),
          State(7, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 1.),
          State(3, Parity::unknown)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, -1.),
          State(7, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, -1.),
          State(3, Parity::unknown)}}),
    new W_pol_dir(
        State(3, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 100.),
          State(7, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 100.),
          State(3, Parity::unknown)}}),

    // 5/2 -> 3/2 -> 5/2
    new W_dir_dir(
        State(5, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(3, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(3, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 1.),
          State(3, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, -1.),
          State(3, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 100.),
          State(3, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(5, Parity::unknown)}}),

    // 5/2 -> 5/2 -> 5/2
    new W_dir_dir(
        State(5, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(5, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 1.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, -1.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 100.),
          State(5, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(5, Parity::unknown)}}),

    // 5/2 -> 7/2 -> 5/2
    new W_dir_dir(
        State(5, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(7, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(7, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 1.),
          State(7, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, -1.),
          State(7, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, -1.),
          State(5, Parity::unknown)}}),
    new W_pol_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 100.),
          State(7, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(5, Parity::unknown)}}),

    // 5/2 -> 9/2 -> 5/2
    new W_dir_dir(
        State(5, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(9, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(5, Parity::unknown)}}),
    new W_dir_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(9, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(5, Parity::unknown)}}),
    new W_dir_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 1.),
          State(9, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 1.),
          State(5, Parity::unknown)}}),
    new W_dir_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, -1.),
          State(9, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, -1.),
          State(5, Parity::unknown)}}),
    new W_dir_dir(
        State(5, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 100.),
          State(9, Parity::positive)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 100.),
          State(5, Parity::unknown)}}),

    // 0 -> 1 -> 1 -> 0
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::unknown)}}),
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::unknown)}}),
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::unknown)}}),

    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::magnetic, 2, EMCharacter::electric, 4, 0.),
          State(2, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(2, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(0, Parity::unknown)}}),

    // 0 -> 2 -> 2 -> 0
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(0, Parity::unknown)}}),
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(0, Parity::unknown)}}),
    new W_dir_dir(
        State(0, Parity::unknown),
        {{Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(0, Parity::unknown)}}),

    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 0.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(0, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 1.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(0, Parity::unknown)}}),
    new W_pol_dir(
        State(0, Parity::positive),
        {{Transition(EMCharacter::electric, 4, EMCharacter::magnetic, 6, 0.),
          State(4, Parity::positive)},
         {Transition(EMCharacter::unknown, 2, EMCharacter::unknown, 4, 100.),
          State(4, Parity::unknown)},
         {Transition(EMCharacter::unknown, 4, EMCharacter::unknown, 6, 0.),
          State(0, Parity::unknown)}}),
};

/**
 * \brief Test upper limits of angular correlations.
 *
 * Test by sampling values from different angular correlations and comparing
 * them to the upper limit. The almost uniform sampling is performed using the
 * SpherePointSampler.
 */
int main() {

  const SpherePointSampler sph_poi_sam;

  const unsigned int n_sphere_points = 1000;
  array<vector<double>, 2> sphere_points = sph_poi_sam.sample(n_sphere_points);

  double ang_cor_max{0.}, ang_cor_val{0.}, ang_cor_upp_lim{0.};

  for (auto ang_cor : ang_corrs) {
    ang_cor_max = 0.;
    ang_cor_upp_lim = ang_cor->get_upper_limit();

    for (size_t i = 0; i < n_sphere_points; ++i) {
      ang_cor_val =
          ang_cor->operator()(sphere_points[0][i], sphere_points[1][i]);
      assert(ang_cor_val <= ang_cor_upp_lim);
      if (ang_cor_val > ang_cor_max) {
        ang_cor_max = ang_cor_val;
      }
    }
  }
}
