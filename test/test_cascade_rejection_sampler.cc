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

    Copyright (C) 2021 Udo Friman-Gayer
*/

#include <array>

using std::array;

#include <vector>

using std::vector;

#include <cassert>

#include <gsl/gsl_math.h>

#include "AngCorrRejectionSampler.hh"
#include "AngularCorrelation.hh"
#include "CascadeRejectionSampler.hh"
#include "EulerAngleRotation.hh"
#include "SphereRejectionSampler.hh"
#include "TestUtilities.hh"

/**
 * Test the cascade rejection sampler.
 * In this test, a cascade of two gamma rays is simulated using the CascadeRejectionSampler.
 * The same cascade is manually reconstructed using two AngCorrRejectionSampler objects and
 * a SphereRejectionSampler to obtain the random vectors and EulerAngleRotation to rotate the
 * random vectors into the expected orientation.
 * This test involves some random-number sampling.
 * By initializing all samplers with the same seed, it is ensured that all calculations give
 * the same results.
 */
int main(){

    const double epsilon = 1e-8;

    // Create a cascade of two angular correlations.
    vector<AngularCorrelation> cascade{
        AngularCorrelation(
            State(0, positive),
            {
                {Transition(magnetic, 2, electric, 4, 0.), State(2, positive)},
                {Transition(magnetic, 2, electric, 4, 0.), State(4, positive)}
            }
        ),
        AngularCorrelation(
            State(2),
            {
                {Transition(2, 4, 0.), State(4)},
                {Transition(4, 6, 0.), State(0)}
            }
        )
    };

    AngCorrRejectionSampler ang_cor_rej_sam_1(cascade[0], 0);
    AngCorrRejectionSampler ang_cor_rej_sam_2(cascade[1], 0);

    CascadeRejectionSampler cas_rej_sam_random(cascade, 0); // Cascade sampler with random initial orientation.
    CascadeRejectionSampler cas_rej_sam_z_axis(cascade, 0, {0., 0., 0.}, true); // Cascade sampler with orientation along the z axis (i.e. no rotation)

    // Test the cascade sampler along the z axis.
    const array<double, 2> theta_phi_single_1 = ang_cor_rej_sam_1(); // Sample first vector from AngCorrRejectionSampler
    const array<double, 2> theta_phi_single_2 = ang_cor_rej_sam_2(); // Sample second vector from AngCorrRejectionSampler
    const vector<array<double, 2>> theta_phi_z_axis = cas_rej_sam_z_axis(); // Sample both from CascadeRejectionSampler

    //      Initial direction of the cascade
    test_numerical_equality<double>(0., theta_phi_z_axis[0][0], epsilon);
    // Since the angle theta is zero, the angle phi is undefined.
    // Internally, the angle phi for the initial direction is inferred from the Euler angle Psi,
    // which gives phi = pi/2 - Psi = pi/2 - 0 = pi/2.
    test_numerical_equality<double>(M_PI_2, theta_phi_z_axis[0][1], epsilon);

    //      First step of the cascade
    test_numerical_equality<double>(theta_phi_single_1[0], theta_phi_z_axis[1][0], epsilon);
    test_numerical_equality<double>(theta_phi_single_1[1], theta_phi_z_axis[1][1], epsilon);

    //      Second step of the cascade
    // This vector must be rotated into the reference frame given by the first emission.
    const EulerAngleRotation eul_ang_rot;
    array<double, 2> theta_phi_single_2_rotated = eul_ang_rot.rotate(theta_phi_single_2, {0., theta_phi_single_1[0], -theta_phi_single_1[1]+M_PI_2});
    test_numerical_equality<double>(theta_phi_single_2_rotated[0], theta_phi_z_axis[2][0], epsilon);
    test_numerical_equality<double>(theta_phi_single_2_rotated[1], theta_phi_z_axis[2][1], epsilon);

    // Test the cascade sampler with a random orientation.
    // Here, both vectors must be rotated into another reference frame.
    SphereRejectionSampler sph_rej_sam([]([[maybe_unused]] const double theta, [[maybe_unused]] const double phi){return 1.;}, 1., 0);
    array<double, 2> sphere_random_vector = sph_rej_sam();
    const vector<array<double, 2>> theta_phi_random = cas_rej_sam_random();
    //      First step of the cascade
    const array<double, 2> theta_phi_single_1_rotated = eul_ang_rot.rotate(theta_phi_single_1, {0., sphere_random_vector[0], -sphere_random_vector[1]+M_PI_2});

    test_numerical_equality<double>(theta_phi_single_1_rotated[0], theta_phi_random[1][0], epsilon);
    test_numerical_equality<double>(theta_phi_single_1_rotated[1], theta_phi_random[1][1], epsilon); 

    //      Second step of the cascade
    theta_phi_single_2_rotated = eul_ang_rot.rotate(theta_phi_single_2, {0., theta_phi_single_1_rotated[0], -theta_phi_single_1_rotated[1]+M_PI_2});

    test_numerical_equality<double>(theta_phi_single_2_rotated[0], theta_phi_random[2][0], epsilon);
    test_numerical_equality<double>(theta_phi_single_2_rotated[1], theta_phi_random[2][1], epsilon);
}