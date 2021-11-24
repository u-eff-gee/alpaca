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

#include <fstream>

using std::ofstream;

#include <iostream>

using std::endl;

#include <sstream>

using std::stringstream;

#include "State.hh"
#include "Transition.hh"
#include "W_dir_dir.hh"
#include "W_pol_dir.hh"

vector<W_gamma_gamma*> w_gamma_gamma{ 
    new W_dir_dir(
		State(0),
		{
			{Transition(2, 4, 0.), State(2)},
			{Transition(2, 4, 0.), State(4)},
		}        
    ),
    new W_dir_dir(
		State(0),
		{
			{Transition(2, 4, 0.), State(2)},
			{Transition(2, 4, 0.), State(2)},
			{Transition(2, 4, 0.), State(4)},
		}        
    ),
    new W_pol_dir(
		State(0, positive),
		{
			{Transition(magnetic, 2, electric, 4, 0.), State(2, positive)},
			{Transition(magnetic, 2, electric, 4, 0.), State(4, positive)},
		}        
    ),
    new W_pol_dir(
		State(3, negative),
		{
			{Transition(magnetic, 2, electric, 4, 0.), State(5, negative)},
			{Transition(magnetic, 2, electric, 4, 0.), State(3, negative)},
		}        
    ),
    new W_pol_dir(
		State(0, positive),
		{
			{Transition(magnetic, 2, electric, 4, 0.), State(2, positive)},
			{Transition(magnetic, 2, electric, 4, 0.), State(2, positive)},
			{Transition(magnetic, 2, electric, 4, 0.), State(4, positive)},
		}        
    ),
};

int main(){

    const unsigned int precision = 8;

    ofstream texfile("test.tex");
    stringstream texfile_buffer;
    texfile_buffer << "\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n";

    for(auto w: w_gamma_gamma){
        texfile_buffer << "\\begin{equation}\n";
        texfile_buffer << w->get_initial_state().str_rep();
        for(auto cascade_step: w->get_cascade_steps()){
            texfile_buffer << " \\rightarrow " << cascade_step.second.str_rep();
        }
        texfile_buffer << "\n\\end{equation}\n";
        texfile_buffer  << "\\begin{align*}\n" 
                        << w->string_representation() 
                        << "\n\\end{align*}\n";
        texfile_buffer  << "\\begin{align*}\n" 
                        << w->string_representation(precision, {"\\theta", "\\varphi", "\\delta_1", "\\delta_2", "\\delta_3"}) 
                        << "\n\\end{align*}\n\\newpage";
    }
    
    texfile_buffer << "\\end{document}\n";
    texfile << texfile_buffer.str();
    texfile.close();

}