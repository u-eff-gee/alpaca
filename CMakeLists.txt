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
#    Copyright (C) 2021-2023 Udo Friman-Gayer

cmake_minimum_required(VERSION 3.16)
project(alpaca VERSION 1.0.5)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

option(BUILD_DOCUMENTATION "Build doxygen documentation for some classes of alpaca" OFF)
if(BUILD_DOCUMENTATION)
        find_package(Doxygen)
        configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in Doxyfile)
        add_custom_target(doc_doxygen ALL COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile)
endif(BUILD_DOCUMENTATION)

option(BUILD_TESTS "Build self tests." OFF)
if(BUILD_TESTS)
        add_compile_options(-fPIC -Wall -Wextra -Wpedantic -fprofile-arcs -ftest-coverage --coverage)
        link_libraries(gcov)
endif(BUILD_TESTS)

add_compile_options(-fPIC -Wall -Wextra -Wpedantic)

include(GNUInstallDirs)

add_subdirectory(python)
add_subdirectory(source)
if(BUILD_TESTS)
        add_subdirectory(test)
endif(BUILD_TESTS)

set(installable_libs angcorrRejectionSampler angular_correlation alphavCoefficient avCoefficient cascadeSampler referenceFrameSampler fCoefficient kappa_coefficient sphereRejectionSampler state stringRepresentable transition uvCoefficient w_dir_dir w_gamma_gamma w_pol_dir)
install(
    TARGETS ${installable_libs}
    EXPORT ALPACA
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT ALPACA
        FILE alpaca-config.cmake
        DESTINATION lib/cmake/alpaca
)
