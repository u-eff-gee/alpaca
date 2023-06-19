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

#include <memory>
#include <vector>

#include <gsl/gsl_math.h>

#include "alpaca/CascadeSampler.hh"
#include "alpaca/DeterministicReferenceFrameSampler.hh"
#include "alpaca/TestUtilities.hh"
#include "alpaca/EulerAngleRotation.hh"

using std::make_shared;
using std::shared_ptr;
using std::vector;

using alpaca::CascadeSampler;
using alpaca::DeterministicReferenceFrameSampler;
using alpaca::ReferenceFrameSampler;
using alpaca::test_numerical_equality;
using alpaca::EulerAngles;

int main() {
  const double epsilon = 1e-8;

  CascadeSampler cascade_sampler(vector<shared_ptr<ReferenceFrameSampler>>{
      make_shared<DeterministicReferenceFrameSampler>(
                                         EulerAngles{0, -M_PI_2, 0}),
      make_shared<DeterministicReferenceFrameSampler>(
                                         EulerAngles{0, -M_PI_2, 0})});

    vector<EulerAngles> cascade = cascade_sampler();

  test_numerical_equality<double>(
            3, cascade[0].data(), EulerAngles{0, -M_PI_2, 0}.data(), epsilon);
  test_numerical_equality<double>(3, cascade[1].data(),
                                  EulerAngles{0, M_PI, 0}.data(), epsilon);
}
