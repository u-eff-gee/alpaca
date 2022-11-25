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

    Copyright (C) 2021, 2022 Udo Friman-Gayer
*/

#include <cassert>

#include "TestUtilities.hh"

int main(){

    bool error_thrown = false;

    test_numerical_equality<double>(1.000, 1.001, 1e-3);

    try{
        test_numerical_equality<double>(1.000, 1.001, 1e-4);
    } catch(const runtime_error &e){
        error_thrown = true;
    }

    assert(error_thrown);
}