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

#include <cmath>

#include "EulerAngleRotation.hh"

array<double, 3> EulerAngleRotation::rotate(const array<double, 3> x_y_z, const array<double, 3> Phi_Theta_Psi) const {

    if(no_rotation_required(Phi_Theta_Psi)){
        return x_y_z;
    }

    const array<array<double, 3>, 3> A = rotation_matrix(Phi_Theta_Psi);

    return array<double, 3>{
        A[0][0]*x_y_z[0] + A[0][1]*x_y_z[1] + A[0][2]*x_y_z[2],
        A[1][0]*x_y_z[0] + A[1][1]*x_y_z[1] + A[1][2]*x_y_z[2],
        A[2][0]*x_y_z[0] + A[2][1]*x_y_z[1] + A[2][2]*x_y_z[2]
    };
}

array<double, 2> EulerAngleRotation::rotate(const array<double, 2> theta_phi, const array<double, 3> Phi_Theta_Psi) const {

    if(no_rotation_required(Phi_Theta_Psi)){
        return theta_phi;
    }

    return get_theta_phi(rotate(get_x_y_z_norm(theta_phi), Phi_Theta_Psi));
}

array<double, 3> EulerAngleRotation::rotate_back(const array<double, 3> xp_yp_zp, const array<double, 3> Phi_Theta_Psi) const {

    if(no_rotation_required(Phi_Theta_Psi)){
        return xp_yp_zp;
    }
 
    return rotate(xp_yp_zp, {-Phi_Theta_Psi[2], -Phi_Theta_Psi[1], -Phi_Theta_Psi[0]});
}

array<double, 2> EulerAngleRotation::rotate_back(const array<double, 2> thetap_phip, const array<double, 3> Phi_Theta_Psi) const {

    if(no_rotation_required(Phi_Theta_Psi)){
        return thetap_phip;
    }

    return get_theta_phi(rotate_back(get_x_y_z_norm(thetap_phip), Phi_Theta_Psi));
}

array<double, 2> EulerAngleRotation::get_theta_phi(const array<double, 3> x_y_z_norm) const {

    // Convention for the special case where both the x and y component of the vector are zero.
    if(x_y_z_norm[0] == 0. && x_y_z_norm[1] == 0.){
        if(x_y_z_norm[2] >= 0.){
            return {0., 0.};
        }
        return {M_PI, 0.};
    }

    const double theta = acos(x_y_z_norm[2]);
    if(x_y_z_norm[0] != 0.){
        if(x_y_z_norm[1] == 0.){
            return {theta, x_y_z_norm[0] > 0. ? 0. : M_PI};
        }
        return {theta, atan(x_y_z_norm[1]/x_y_z_norm[0])};
    }

    return {theta, atan(x_y_z_norm[1] >=0 ? INFINITY : -INFINITY)};
}

array<double, 3> EulerAngleRotation::get_x_y_z_norm(const array<double, 2> theta_phi) const {
    
    double cos_the{cos(theta_phi[0])}, sin_the{sin(theta_phi[0])}, cos_phi{cos(theta_phi[1])}, sin_phi{sin(theta_phi[1])};

    return array<double, 3>{ 
        sin_the*cos_phi,
        sin_the*sin_phi,
        cos_the
    };

}

array<array<double, 3>, 3> EulerAngleRotation::rotation_matrix(const array<double, 3> Phi_Theta_Psi) const {
    
    const double cos_phi{cos(Phi_Theta_Psi[0])}, sin_phi{sin(Phi_Theta_Psi[0])}, cos_the{cos(Phi_Theta_Psi[1])}, sin_the{sin(Phi_Theta_Psi[1])}, cos_psi{cos(Phi_Theta_Psi[2])}, sin_psi{sin(Phi_Theta_Psi[2])};

    return array<array<double, 3>, 3> {
        array<double, 3>{ 
            cos_psi*cos_phi - cos_the*sin_phi*sin_psi,
            cos_psi*sin_phi + cos_the*cos_phi*sin_psi,
            sin_psi*sin_the
        },
        array<double, 3>{
            -sin_psi*cos_phi - cos_the*sin_phi*cos_psi,
            -sin_psi*sin_phi + cos_the*cos_phi*cos_psi,
            cos_psi*sin_the
        },
        array<double, 3>{
            sin_the*sin_phi,
            -sin_the*cos_phi,
            cos_the
        }
    };

}

bool EulerAngleRotation::no_rotation_required(const array<double, 3> Phi_Theta_Psi) const {
    if(Phi_Theta_Psi[0] == 0. && Phi_Theta_Psi[1] == 0. && Phi_Theta_Psi[2] == 0.){
        return true;
    }
    return false;
}