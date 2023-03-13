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

using std::array;

#include "AngCorrRejectionSampler.hh"
#include "EulerAngleRotation.hh"

AngCorrRejectionSampler::AngCorrRejectionSampler(AngularCorrelation &w,
                                                 const int seed,
                                                 const unsigned int max_tri)
    : SphereRejectionSampler(w, w.get_upper_limit(), seed, max_tri) {}