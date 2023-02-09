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

#include "AngularCorrelation.hh"
#include "State.hh"
#include "Transition.hh"

void calculate_angular_correlation(){
    AngularCorrelation ang_corr(
        State(0, positive), 
        {
            {
                Transition(electric, 2, magnetic, 4, 0.),
                State(2, negative)
            },
            {
                Transition(electric, 2, magnetic, 4, 0.),
                State(0, positive)
            }
        }
    );
}

/**
 * This code can be used to test for memory leaks of alpaca.
 * A common use case of alpaca is the repeated construction and deletion of AngularCorrelation 
 * objects, which is particularly important for the python interface.
 * 
 * The ang_corr object created in the function above should be destroyed completely when it goes 
 * out of scope at the end of the function call, i.e. after some short initialization, this program
 * should not allocate additional memory any more.
 * 
 * By profiling this test with valgrind, the author was able to find a memory leak in 
 * AngularCorrelation.
 * It was decided to simply keep the script as an additional test.
 */
int main(){
    for(int i = 0; i < 10000; ++i){
        calculate_angular_correlation();
    }
}