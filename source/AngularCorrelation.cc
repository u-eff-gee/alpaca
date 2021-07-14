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

#include <stdexcept>
#include <string>

#include "AngularCorrelation.hh"
#include "TestUtilities.hh"
#include "W_dir_dir.hh"
#include "W_pol_dir.hh"

using std::invalid_argument;
using std::to_string;

AngularCorrelation::AngularCorrelation(const State ini_sta, const vector<pair<Transition, State>> cas_ste):
    euler_angle_rotation(EulerAngleRotation()), w_gamma_gamma(nullptr)
{
    check_cascade(ini_sta, cas_ste);

    if(cas_ste[0].first.em_char == em_unknown){
        w_gamma_gamma = std::make_unique<W_dir_dir>(ini_sta, cas_ste);
    } else{
        w_gamma_gamma = std::make_unique<W_pol_dir>(ini_sta, cas_ste);
    }
}

AngularCorrelation::AngularCorrelation(const State ini_sta, const vector<State> cas_sta):
    euler_angle_rotation(EulerAngleRotation()), w_gamma_gamma(nullptr){

    if(cas_sta.size() < 2){
        throw invalid_argument("Cascade must have at least two transition - state pairs.");
    }

    vector<pair<Transition, State>> cascade_steps;

    cascade_steps.push_back(
        {
            infer_transition({ini_sta, cas_sta[0]}),
            cas_sta[0]
        }
    );

    for(size_t i = 0; i < cas_sta.size()-1; ++i){
        cascade_steps.push_back(
            {
                infer_transition({cas_sta[i], cas_sta[i+1]}),
                cas_sta[i+1]
            }
        );
    }

    if(cascade_steps[0].first.em_char == em_unknown){
        w_gamma_gamma = std::make_unique<W_dir_dir>(ini_sta, cascade_steps);
    } else{
        w_gamma_gamma = std::make_unique<W_pol_dir>(ini_sta, cascade_steps);
    }

}

Transition AngularCorrelation::infer_transition(const pair<State, State> states) const {
    
    if(states.first.two_J == 0 and states.second.two_J == 0){
        throw invalid_argument("An electromagnetic transition between two spin-0 states with the absorption/emission of a single photon is not possible.");
    }

    int two_L = abs(states.first.two_J - states.second.two_J);
    if(two_L == 0){
        two_L += 2;
    }
    EMCharacter em = em_unknown;
    EMCharacter emp = em_unknown;
    if(states.first.parity != parity_unknown && states.second.parity != parity_unknown){
        if(two_L % 4 == 0){
            if(states.first.parity == states.second.parity){
                em = electric;
                emp = magnetic;
            } else{
                em = magnetic;
                emp = electric;
            }
        } else{
            if(states.first.parity == states.second.parity){
                em = magnetic;
                emp = electric;
            } else{
                em = electric;
                emp = magnetic;
            }
        }
    }

    return Transition(em, two_L, emp, two_L + 2, 0.);
}

double AngularCorrelation::operator()(const double theta, const double phi) const {
    return w_gamma_gamma->operator()(theta, phi);
}

double AngularCorrelation::operator()(const double theta, const double phi, const array<double, 3> euler_angles) const {
    array<double, 2> thetap_phip = euler_angle_rotation.rotate_back(array<double, 2>{theta, phi}, euler_angles);

    return (*this)(thetap_phip[0], thetap_phip[1]);
}

void AngularCorrelation::check_cascade(const State ini_sta, const vector<pair<Transition, State>> cas_ste) const {

    if(cas_ste.size() < 2){
        throw invalid_argument("Cascade must have at least two transition - state pairs.");
    }

    check_angular_momenta(ini_sta, cas_ste);

    check_triangle_inequalities(ini_sta, cas_ste);

    check_em_transitions(ini_sta, cas_ste);

}

void AngularCorrelation::check_angular_momenta(const State ini_sta, const vector<pair<Transition, State>> cas_ste) const {

    const int even_odd = ini_sta.two_J % 2;

    for(size_t i = 0; i < cas_ste.size(); ++i){
        if(cas_ste[i].second.two_J % 2 != even_odd){
            throw invalid_argument("Unphysical mixing of half-integer and integer spins in cascade.");
        }
    }
}

void AngularCorrelation::check_triangle_inequalities(const State ini_sta, const vector<pair<Transition, State>> cas_ste) const {
    
    if(
        !fulfils_triangle_inequality<int>(ini_sta.two_J, cas_ste[0].second.two_J, cas_ste[0].first.two_L) && !fulfils_triangle_inequality<int>(ini_sta.two_J, cas_ste[0].second.two_J, cas_ste[0].first.two_Lp)
    ){
        throw invalid_argument("Triangle inequality selection rule not fulfilled for any multipolarity of transition #1: " + cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
    }

    for(size_t i = 1; i < cas_ste.size(); ++i){
        if(
            !fulfils_triangle_inequality<int>(cas_ste[i-1].second.two_J, cas_ste[i].second.two_J, cas_ste[i].first.two_L) && !fulfils_triangle_inequality<int>(cas_ste[i-1].second.two_J, cas_ste[i].second.two_J, cas_ste[i].first.two_Lp)
        ){
            throw invalid_argument("Triangle inequality selection rule not fulfilled for any multipolarity of transition #" + cas_ste[i].first.str_rep(cas_ste[i-1].second, cas_ste[i].second));
        }
    }
}

void AngularCorrelation::check_em_transitions(const State ini_sta, const vector<pair<Transition, State>> cas_ste) const {

    if(ini_sta.parity != parity_unknown && cas_ste[0].second.parity != parity_unknown){
        if(cas_ste[0].first.em_char != em_unknown){
            if(!valid_em_character(ini_sta.parity, cas_ste[0].second.parity, cas_ste[0].first.two_L, cas_ste[0].first.em_char)){
                throw invalid_argument("Incorrect electromagnetic character '" + Transition::em_str_rep(cas_ste[0].first.em_char) + "' for transition #1: " + cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
            }

            if(cas_ste[0].first.em_charp != em_unknown){
                if(!valid_em_character(ini_sta.parity, cas_ste[0].second.parity, cas_ste[0].first.two_Lp, cas_ste[0].first.em_charp)){
                    throw invalid_argument("Incorrect electromagnetic character '" + Transition::em_str_rep(cas_ste[0].first.em_charp) + "' for transition #1: " + cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
                }
            } else{
                throw invalid_argument("Only one electromagnetic character defined for transition #1: " + cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));
            }
        }

        if(cas_ste[0].first.em_char == em_unknown && cas_ste[0].first.em_charp != em_unknown){
            throw invalid_argument("Only one electromagnetic character defined for transition #1: " + cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));        
        }
    } else if(cas_ste[0].first.em_char != em_unknown || cas_ste[0].first.em_charp != em_unknown){
        throw invalid_argument("Electromagnetic character defined, but one or both parities missing for transition #1: " + cas_ste[0].first.str_rep(ini_sta, cas_ste[0].second));        
    }

    for(size_t i = 1; i < cas_ste.size(); ++i){
        if(cas_ste[i-1].second.parity != parity_unknown && cas_ste[i].second.parity != parity_unknown){
            if(cas_ste[i].first.em_char != em_unknown){    
                if(!valid_em_character(cas_ste[i-1].second.parity, cas_ste[i].second.parity, cas_ste[i].first.two_L, cas_ste[i].first.em_char)){
                    throw invalid_argument("Incorrect electromagnetic character '" + Transition::em_str_rep(cas_ste[0].first.em_char) + "' for transition #" + to_string(i+1) + ": " + cas_ste[i].first.str_rep(cas_ste[i-1].second, cas_ste[i].second));

                }
            }

            if(cas_ste[i].first.em_charp != em_unknown){    
                if(!valid_em_character(cas_ste[i-1].second.parity, cas_ste[i].second.parity, cas_ste[i].first.two_Lp, cas_ste[i].first.em_charp)){
                    throw invalid_argument("Incorrect electromagnetic character '" + Transition::em_str_rep(cas_ste[0].first.em_charp) + "' for transition #" + to_string(i+1) + ": " + cas_ste[i].first.str_rep(cas_ste[i-1].second, cas_ste[i].second));

                }
            } else{
                throw invalid_argument("Only one electromagnetic character defined for transition #"
                + to_string(i+1) + ": " + cas_ste[i].first.str_rep(cas_ste[i-1].second, cas_ste[i].second)
                );
            }

            if(cas_ste[i].first.em_char == em_unknown && cas_ste[i].first.em_charp != em_unknown){
                throw invalid_argument("Only one electromagnetic character defined for transition #" + to_string(i+1) + ": " + cas_ste[i].first.str_rep(cas_ste[i-1].second, cas_ste[i].second));        
            }
        } else if(cas_ste[i].first.em_char != em_unknown || cas_ste[i].first.em_charp != em_unknown){
            throw invalid_argument("Electromagnetic character defined, but one or both parities missing for transition #" + to_string(i+1) + ": " + cas_ste[i].first.str_rep(cas_ste[i-1].second, cas_ste[i].second));
        }        
    }
}

bool AngularCorrelation::valid_em_character(const Parity p0, const Parity p1, const int two_L, const EMCharacter em) const {

    if(p0 == p1){
        if((two_L/2) % 2 == 0){
            if(em != electric){
                return false;
            }
        } else if(em != magnetic){
            return false;
        }

        return true;
    }

    if((two_L/2) % 2 == 0){
        if(em != magnetic){
            return false;
        }
    } else if(em != electric){
        return false;
    }

    return true;

}

extern "C" {
    double angular_correlation(const double theta, const double phi, const size_t n_cas_ste, int* two_J, short* par, short* em_char, int* two_L, short* em_charp, int* two_Lp, double* delta, double* PhiThetaPsi){
        State initial_state{two_J[0], (Parity) par[0]};
        vector<pair<Transition, State>> cascade_steps;

        for(size_t i = 0; i < n_cas_ste; ++i){
            cascade_steps.push_back(
                {
                    Transition{(EMCharacter) em_char[i], two_L[i], (EMCharacter) em_charp[i], two_Lp[i], delta[i]},
                    State{two_J[i+1], (Parity) par[i+1]}
                }
            );
        }

        return AngularCorrelation(initial_state, cascade_steps).operator()(theta, phi, {PhiThetaPsi[0], PhiThetaPsi[1], PhiThetaPsi[2]});
    }

    void* create_angular_correlation(const size_t n_cas_ste, int* two_J, short* par, short* em_char, int* two_L, short* em_charp, int* two_Lp, double* delta){

        State initial_state{two_J[0], (Parity) par[0]};
        vector<pair<Transition, State>> cascade_steps;

        for(size_t i = 0; i < n_cas_ste; ++i){
            cascade_steps.push_back(
                {
                    Transition{(EMCharacter) em_char[i], two_L[i], (EMCharacter) em_charp[i], two_Lp[i], delta[i]},
                    State{two_J[i+1], (Parity) par[i+1]}
                }
            );
        }

        return new AngularCorrelation(initial_state, cascade_steps);
    }

    void* create_angular_correlation_with_transition_inference(const size_t n_cas_ste, int* two_J, short* par){

        State initial_state{two_J[0], (Parity) par[0]};
        vector<State> cascade_states;

        for(size_t i = 0; i < n_cas_ste; ++i){
            cascade_states.push_back(
                {
                    State{two_J[i+1], (Parity) par[i+1]}
                }
            );
        }

        return new AngularCorrelation(initial_state, cascade_states);
    }

    void evaluate_angular_correlation(AngularCorrelation *angular_correlation, const size_t n_angles, double* theta, double* phi, double* result){

        for(size_t i = 0; i < n_angles; ++i){
            result[i] = angular_correlation->operator()(theta[i], phi[i]);
        }

    }

    void evaluate_angular_correlation_rotated(AngularCorrelation *angular_correlation, const size_t n_angles, double* theta, double* phi, double* PhiThetaPsi, double* result){

        for(size_t i = 0; i < n_angles; ++i){
            result[i] = angular_correlation->operator()(theta[i], phi[i], {PhiThetaPsi[0], PhiThetaPsi[1], PhiThetaPsi[2]});
        }

    }

    void free_angular_correlation(AngularCorrelation *angular_correlation){
        delete angular_correlation;
    }

    void get_em_char(AngularCorrelation *angular_correlation, short* em_char){
        
        vector<pair<Transition, State>> cascade_steps = angular_correlation->get_cascade_steps();
        for(size_t i = 0; i < cascade_steps.size(); ++i){
            em_char[i] = cascade_steps[i].first.em_char;
        }

    }

    void get_two_L(AngularCorrelation *angular_correlation, int* two_L){
        
        vector<pair<Transition, State>> cascade_steps = angular_correlation->get_cascade_steps();
        for(size_t i = 0; i < cascade_steps.size(); ++i){
            two_L[i] = cascade_steps[i].first.two_L;
        }

    }
}