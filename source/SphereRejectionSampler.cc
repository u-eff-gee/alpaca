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

#include <numeric>
#include <utility>
#include <vector>

using std::accumulate;
using std::get;
using std::vector;

#include <gsl/gsl_math.h>

#include "SphereRejectionSampler.hh"

SphereRejectionSampler::SphereRejectionSampler(function<double(const double, const double)> dis, const double dis_max, const int seed, const unsigned int max_tri):
    distribution(dis),
    distribution_maximum(dis_max),
    max_tries(max_tri),
    euler_angle_rotation(EulerAngleRotation())
{
    random_engine = mt19937(seed);
}

pair<unsigned int, array<double, 2>> SphereRejectionSampler::sample(){
    
    array<double, 2> theta_phi;
    double dis_val;

    for(unsigned int i = 0; i < max_tries; ++i){

        theta_phi = sample_theta_phi();
        dis_val = uniform_random(random_engine)*distribution_maximum;

        if(dis_val <= distribution(theta_phi[0], theta_phi[1])){
            return {i+1, {theta_phi[0], theta_phi[1]}};
        }

    }

    return {max_tries, {0., 0.}};
}

array<double, 2> SphereRejectionSampler::operator()(){
    pair<unsigned int, array<double, 2>> sampled_theta_phi = sample();

    return {sampled_theta_phi.second[0], sampled_theta_phi.second[1]};
}

double SphereRejectionSampler::estimate_efficiency(const unsigned int n_tries){
    vector<unsigned int> required_tries(n_tries);
    
    pair<unsigned int, array<double, 2>> sampled_theta_phi;

    for(unsigned int i = 0; i < n_tries; ++i){
        sampled_theta_phi = sample();
        required_tries[i] = sampled_theta_phi.first;
    }

    return (double) n_tries / (double) accumulate(required_tries.begin(), required_tries.end(), 0);
}

double SphereRejectionSampler::sample_theta(){
    return acos(2.*uniform_random(random_engine)-1.);
}

double SphereRejectionSampler::sample_phi(){
    return 2.*M_PI*uniform_random(random_engine);
}

array<double, 2> SphereRejectionSampler::sample_theta_phi(){
    return {sample_theta(), sample_phi()};
}