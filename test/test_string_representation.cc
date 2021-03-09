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

vector<W_dir_dir> w_dir_dir{
    {
		State(0), 
		{
			{Transition(2, 4, 0.), State(2)},
			{Transition(0, 4, 0.), State(0)},
		} 
    },
    {
		State(0), 
		{
			{Transition(2, 4, 0.), State(2)},
			{Transition(2, 4, 0.), State(2)},
		} 
    },
    {
		State(0), 
		{
			{Transition(2, 4, 0.), State(2)},
			{Transition(2, 4, 0.), State(4)},
		} 
    },
    {
		State(0), 
		{
			{Transition(2, 4, 0.), State(2)},
			{Transition(4, 6, 0.), State(6)},
		} 
    },
    {
		State(0), 
		{
			{Transition(4, 6, 0.), State(4)},
			{Transition(4, 6, 0.), State(0)},
		} 
    },
    {
		State(0), 
		{
			{Transition(4, 6, 0.), State(4)},
			{Transition(2, 4, 0.), State(4)},
		} 
    },
    {
		State(0), 
		{
			{Transition(4, 6, 0.), State(4)},
			{Transition(2, 4, 0.), State(6)},
		} 
    },
    {
		State(0), 
		{
			{Transition(4, 6, 0.), State(4)},
			{Transition(4, 6, 0.), State(8)},
		} 
    },
};

vector<W_pol_dir> w_pol_dir = {
    {
		State(0, positive), 
		{
			{Transition(magnetic, 2, electric, 4, 0.), State(2, positive)},
			{Transition(magnetic, 0, electric, 4, 0.), State(0, positive)},
		} 
    },
    {
		State(0, positive), 
		{
			{Transition(electric, 2, magnetic, 4, 0.), State(2, negative)},
			{Transition(electric, 0, magnetic, 4, 0.), State(0, positive)},
		} 
    },
    {
		State(0, positive), 
		{
			{Transition(magnetic, 2, electric, 4, 0.), State(2, positive)},
			{Transition(magnetic, 0, electric, 4, 0.), State(2, positive)},
		} 
    },
    {
		State(0, positive), 
		{
			{Transition(electric, 2, magnetic, 4, 0.), State(2, negative)},
			{Transition(electric, 0, magnetic, 4, 0.), State(2, positive)},
		} 
    },
};

int main(){

    const unsigned int precision = 8;

    ofstream texfile("test.tex");
    stringstream texfile_buffer;
    texfile_buffer << "\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n";

    for(auto w: w_dir_dir){
        texfile_buffer << "\\begin{equation}\n";
        texfile_buffer << w.get_initial_state().str_rep() << " \\rightarrow " << w.get_cascade_steps()[0].second.str_rep() << " \\rightarrow " << w.get_cascade_steps()[1].second.str_rep();
        texfile_buffer << "\n\\end{equation}\n";
        texfile_buffer << "\\begin{multline*}\n" << w.
        string_representation() << "\n\\end{multline*}\n";
        texfile_buffer << "\\begin{multline*}\n" << w.string_representation(precision, {"\\theta", "\\varphi", "\\delta_1", "\\delta_2"}) << "\n\\end{multline*}\n\\newpage";
    }

    for(auto w: w_pol_dir){
        texfile_buffer << "\\begin{equation}\n";
        texfile_buffer << w.get_initial_state().str_rep() << " \\rightarrow " << w.get_cascade_steps()[0].second.str_rep() << " \\rightarrow " << w.get_cascade_steps()[1].second.str_rep();
        texfile_buffer << "\n\\end{equation}\n";
        texfile_buffer << "\\begin{multline*}\n" << w.
        string_representation() << "\n\\end{multline*}\n";
        texfile_buffer << "\\begin{multline*}\n" << w.string_representation(precision, {"\\theta", "\\varphi", "\\delta_1", "\\delta_2"}) << "\n\\end{multline*}\n\\newpage";
    }
    
    texfile_buffer << "\\end{document}\n";
    texfile << texfile_buffer.str();
    texfile.close();

}