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

#include <cassert>

#include "Transition.hh"

int main(){
    // Test IO of the Transition class.

    bool error_thrown = false;

    // Error: Both multipolarities are the same.
    // Test for both possible constructors.
    try{
        Transition transition(2, 2, 0.);
    } catch (const invalid_argument &e){
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    try{
        Transition transition(electric, 2, magnetic, 2, 0.);
    } catch (const invalid_argument &e){
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Check multipolarity IO

    // Error: Multipolarity smaller than zero.
    try{
        Transition transition(electric, -2, magnetic, 2, 0.);
    } catch (const invalid_argument &e){
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Error: Multipolarity zero.
    try{
        Transition transition(electric, 0, magnetic, 2, 0.);
    } catch (const invalid_argument &e){
        error_thrown = true;
    }

    assert(error_thrown);
    error_thrown = false;

    // Check string representation

    Transition transition(2, 4, 0.);
    try{
        transition.em_str_rep(em_unknown);
    } catch(const runtime_error &e){
        error_thrown = true;
    }

    assert(error_thrown);
}