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
project(
        alpaca
        VERSION 1.0.5
        HOMEPAGE_URL "https://github.com/u-eff-gee/alpaca"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

option(BUILD_DOCUMENTATION "Build doxygen documentation for some classes of alpaca" OFF)
option(BUILD_TESTS "Build self tests" OFF)
option(BUILD_PYTHON "Enable building python bindings" ON)
option(ENABLE_COVERAGE "Enable coverage reporting" OFF)
option(WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

if(BUILD_DOCUMENTATION)
        find_package(Doxygen)
        configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in Doxyfile)
        add_custom_target(doc_doxygen ALL COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile)
endif(BUILD_DOCUMENTATION)

add_library(alpaca_options INTERFACE)
target_compile_features(alpaca_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})
add_library(alpaca::options ALIAS alpaca_options)

target_compile_options(
        alpaca_options
        INTERFACE
        -Wall
        -Wextra
        -Wshadow # warn if a variable declaration shadows one from a parent context
        -Wconversion # warn on type conversions that may lose data
        -Wpedantic # warn if non-standard C++ is used
        -Wnon-virtual-dtor # warn if class with virtual functions has non-virtual destructor
        -Wold-style-cast # warn for c-style casts
        -Wcast-align # warn for potential performance problem casts
        -Wunused # warn on anything being unused
        -Woverloaded-virtual # warn if you overload (not override) a virtual function
        -Wsign-conversion # warn on sign conversions
        -Wnull-dereference # warn if a null dereference is detected
        -Wdouble-promotion # warn if float is implicit promoted to double
        -Wformat=2 # warn on security issues around functions that format output (ie printf)
        -Wmisleading-indentation # warn if indentation implies blocks where blocks do not exist
        -Wduplicated-cond # warn if if / else chain has duplicated conditions
        -Wduplicated-branches # warn if if / else branches have duplicated code
        -Wlogical-op # warn about logical operations being used where bitwise were probably wanted
        -Wuseless-cast # warn if you perform a cast to the same type
        -Wpedantic)

if(ENABLE_COVERAGE)
        target_compile_options(
                alpaca_options
                INTERFACE --coverage -O0 -g -fprofile-arcs -ftest-coverage)
        target_link_libraries(alpaca_options INTERFACE --coverage)
        link_libraries(gcov)
endif()

if(WARNINGS_AS_ERRORS)
        target_compile_options(alpaca_options INTERFACE -Werror)
endif()

include(GNUInstallDirs)


add_subdirectory(source)
if(BUILD_TESTS)
        add_subdirectory(test)
endif(BUILD_TESTS)

set(installable_libs
        coeff
        coeff_static
        angcorr
        angcorr_static
        sampler
        sampler_static
        alpaca_options)
install(
       TARGETS ${installable_libs}
       EXPORT ALPACA
       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
       ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
       PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
       EXPORT ALPACA
       FILE alpaca-config.cmake
       NAMESPACE alpaca::
       DESTINATION lib/cmake/alpaca
)
