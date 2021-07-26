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

#include <cassert>

#include <gsl/gsl_sf.h>

#include "AngCorrRejectionSampler.hh"
#include "AngularCorrelation.hh"
#include "SphereRejectionSampler.hh"
#include "TestUtilities.hh"

/**
 * Test the AngCorrRejectionSampler by verifying that the class, using an AngularCorrelation 
 * object, does exactly the same as a SphereRejectionSampler, using an equivalent analytical
 * expression of the angular correlation.
 */
int main(){

    const int seed = 0;

    AngularCorrelation ang_cor(
        State(0, positive),
        {
            {Transition(magnetic, 2, electric, 4, 0.), State(2, positive)},
            {Transition(magnetic, 2, electric, 4, 0.), State(0, positive)}
        }
    );

    AngCorrRejectionSampler ang_cor_sam(ang_cor, seed);

    // Analytical expression for the angular correlation above.
    // See, e.g., Eq. (1) in Ref. \cite Pietralla2001.
    SphereRejectionSampler sph_rej_sam(
        [](const double theta, const double phi){ return 1.+0.5*(
		gsl_sf_legendre_Pl(2, cos(theta))
		+0.5*cos(2.*phi)*gsl_sf_legendre_Plm(2, 2, cos(theta))); }, ang_cor.get_upper_limit(), seed
    );

    array<double, 2> theta_phi_1;
    array<double, 2> theta_phi_2;

    for(unsigned int n = 0; n < 1000; ++n){
        theta_phi_1 = ang_cor_sam();
        theta_phi_2 = sph_rej_sam();

        test_numerical_equality<double>(theta_phi_1[0], theta_phi_2[0], 1e-6);
        test_numerical_equality<double>(theta_phi_1[1], theta_phi_2[1], 1e-6);
    }

    // Check that the default values 
    //
    // theta = 0
    // phi   = 0
    //
    // are returned when AngCorrRejectionSampler can not find a valid vector.
    // In order to test it, the max_tri is set to zero.
    // This way, the random sampling is bypassed.
    AngCorrRejectionSampler ang_cor_sam_2(ang_cor, 0, 0);
    theta_phi_1 = ang_cor_sam_2();    
    assert(theta_phi_1[0] == 0.);
    assert(theta_phi_1[1] == 0.);
}