#    This file is part of alpaca.
#
#    alpaca is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    alpaca is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with alpaca.  If not, see <https://www.gnu.org/licenses/>.
#
#    Copyright (C) 2021 Udo Friman-Gayer

cmake_minimum_required(VERSION 3.16)
project(alpaca VERSION 0.1.0)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

option(BUILD_DOCUMENTATION "Build doxygen documentation for some classes of alpaca" OFF)
if(BUILD_DOCUMENTATION)
        find_package(Doxygen)
        configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in Doxyfile)
        add_custom_target(doc_doxygen ALL COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile)
endif(BUILD_DOCUMENTATION)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(source)
add_subdirectory(test)