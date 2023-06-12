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

#include <array>

using std::array;

#include "alpaca/ReferenceFrameSampler.hh"

namespace alpaca {

class DeterministicReferenceFrameSampler : public ReferenceFrameSampler {
public:
  DeterministicReferenceFrameSampler(const array<double, 3> Phi_Theta_Psi)
      : Phi_Theta_Psi(Phi_Theta_Psi) {}

  pair<unsigned int, array<double, 3>> sample() override {
    return {1, Phi_Theta_Psi};
  }

protected:
  const array<double, 3> Phi_Theta_Psi;
};

} // namespace alpaca
