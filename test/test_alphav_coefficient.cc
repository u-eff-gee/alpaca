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

#include "AlphavCoefficient.hh"

int main(){
    // As a test case, use the nu=2 coefficient for the polarized 0->1 transition in the 
    // polarization-direction correlation in Sec. "4 Numerical example" of Ref. \cite Iliadis2021.
    // Here, the difference between the Ev coefficients used in the aforementioned publication and 
    // the kappa coefficients used here is once more apparent:
    // In this code, the negative sign due to the magnetic character of the transition is not 
    // included in W_LP but applied later as a factor.
    // Therefore, the entire term seems to have a different sign as in Ref. \cite Iliadis2021.
    // See also the discussion in `alpaca/test/test_Ev_coefficient.cc`.
    AlphavCoefficient av_coef(4, 2, 4, 0, 2);
    const string str_rep = 
        string("(-1)") 
        + KappaCoefficient(4, 2, 2).string_representation() + FCoefficient(4, 2, 2, 0, 2).string_representation()
        + "+2" + KappaCoefficient(4, 2, 4).string_representation() + FCoefficient(4, 2, 4, 0, 2).string_representation()
        + "\\delta"
        + "+" + KappaCoefficient(4, 4, 4).string_representation() + FCoefficient(4, 4, 4, 0, 2).string_representation()
        + "\\delta^{2}";
    assert(av_coef.string_representation() == str_rep);
    assert(av_coef.string_representation(3) == "(-1)\\times\\left(-0.5\\right)\\times0.707+2\\times\\left(-0.167\\right)\\times0\\times\\delta+0.5\\times0\\times\\delta^{2}");

}