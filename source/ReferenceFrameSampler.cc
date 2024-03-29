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

#include <numeric>

using std::accumulate;

#include <utility>

using std::pair;

#include <vector>

using std::vector;

#include "ReferenceFrameSampler.hh"

#include "EulerAngleRotation.hh"

array<double, 3> ReferenceFrameSampler::operator()() { return sample().second; }

double ReferenceFrameSampler::estimate_efficiency(const unsigned int n_tries) {
  vector<unsigned int> required_tries(n_tries);

  for (unsigned int i = 0; i < n_tries; ++i) {
    required_tries[i] = sample().first;
  }

  return (double)n_tries /
         (double)accumulate(required_tries.begin(), required_tries.end(), 0);
}