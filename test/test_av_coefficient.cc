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

#include "AvCoefficient.hh"

int main(){
    // As a test case, use the nu=0 coefficient for the mixed 1->2 transition in the direction-direction correlation in Sec. "4 Numerical example" of Ref. \cite Iliadis2021.
    AvCoefficient av_coef(0, 2, 4, 4, 2);
    const string str_rep = 
        FCoefficient(0, 2, 2, 4, 2).string_representation()
        +"+2"+FCoefficient(0, 2, 4, 4, 2).string_representation() + "\\delta"
        +"+"+FCoefficient(0, 4, 4, 4, 2).string_representation() + "\\delta^{2}";
    assert(av_coef.string_representation() == str_rep);
    assert(av_coef.string_representation(3) == "1+2\\times0\\times\\delta+1\\times\\delta^{2}");

}
