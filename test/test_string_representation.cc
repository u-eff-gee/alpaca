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

#include <iostream>

using std::cout;
using std::endl;

#include "State.hh"
#include "Transition.hh"
#include "W_dir_dir.hh"

int main(){

	W_dir_dir w_dir_dir(
		State(0), 
		{
			{Transition(2, 4, 0.), State(2)},
			{Transition(2, 4, 0.), State(2)},
		} 
	);

    cout << w_dir_dir.string_representation() << endl;

}