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

#pragma once

#include <string>

using std::string;

#include <vector>

using std::vector;

/**'
 * \brief Abstract class for string-representable expressions
 */
class StringRepresentable{

public:
    /**
     * \brief Return string representation of expression.
     * 
     * \param variable_names Names for the variables of a function (default: {} i.e. use default names).
     * 
     * \return String representation.
     */
    virtual string string_representation(const vector<string> variable_names = {}) const = 0;
};