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

#include <cassert>

#include <gsl/gsl_math.h>

#include "EulerAngleRotation.hh"
#include "SphereRejectionSampler.hh"
#include "TestUtilities.hh"

/**
 * This script tests the correctness of SphereRejectionSampler by sampling from a distribution
 * that is 1 for \f$\varphi < \pi\f$ and 0 otherwise.
 * For an optimum choice \f$W_\mathrm{max} = 1\f$, this should result in an efficiency of 
 * \f$\epsilon = 0.5\f$.
 * For the choice \f$W_\mathrm{max} = 2\f$, this should result in an efficiency of 
 * \f$\epsilon = 0.25\f$.
 * 
 * It also tests the Euler-angle rotation functionality of SphereRejectionSampler by sampling two
 * vectors, one from an unrotated distribution and another one from a rotated distribution.
 * By using the same random number seed and two different SphereRejectionSampler objects, it is 
 * ensured that the random vector before the rotation was the same as the unrotated one.
 * By manually applying the rotation using the EulerAngleRotation class, the unrotated vector is
 * transformed into the rotated vector's reference frame and both are tested for equality.
 */
int main(){

    SphereRejectionSampler sph_rej_sam([]([[maybe_unused]] const double theta, const double phi){ return phi < M_PI ? 1. : 0.; }, 1., 0);

    double efficiency = sph_rej_sam.estimate_efficiency(1e5);

    test_numerical_equality<double>(efficiency, 0.5, 1e-3);

    SphereRejectionSampler sph_rej_sam_2([]([[maybe_unused]] const double theta, const double phi){ return phi < M_PI ? 1. : 0.; }, 2., 0);

    efficiency = sph_rej_sam_2.estimate_efficiency(1e5);

    test_numerical_equality<double>(efficiency, 0.25, 1e-3);

    // Test the case in which no vector can be found.
    SphereRejectionSampler sph_rej_sam_3([]([[maybe_unused]] const double theta, [[maybe_unused]] const double phi){ return -1.; }, 0.5, 0);
    const pair<unsigned int, array<double, 2>> theta_phi_default = sph_rej_sam_3.sample();
    assert(theta_phi_default.first == 1000);
    assert(theta_phi_default.second[0] == 0.);
    assert(theta_phi_default.second[1] == 0.);

    // Test Euler-angle rotation

    SphereRejectionSampler sph_rej_sam_uni([]([[maybe_unused]] const double theta, [[maybe_unused]] const double phi){return 1.;}, 1., 0);
    const array<double, 2> theta_phi = sph_rej_sam_uni();
    SphereRejectionSampler sph_rej_sam_uni_2([]([[maybe_unused]] const double theta, [[maybe_unused]] const double phi){return 1.;}, 1., 0);
    const array<double, 3> euler_angles = {0.1, 0.2, 0.3};
    array<double, 2> theta_phi_rotated = sph_rej_sam_uni_2(euler_angles);

    const EulerAngleRotation euler_angle_rotation;
    array<double, 2> theta_phi_rotated_manually = euler_angle_rotation.rotate(theta_phi, euler_angles);

    // Test that the rotation has an effect.
    assert(theta_phi[0] != theta_phi_rotated[0]);
    assert(theta_phi[1] != theta_phi_rotated[1]);
    // Test that the rotation can be reproduced.
    assert(theta_phi_rotated[0] == theta_phi_rotated_manually[0]);
    assert(theta_phi_rotated[1] == theta_phi_rotated_manually[1]);
}