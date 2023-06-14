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

set(libcoeff_src
        AlphavCoefficient.cc
        AvCoefficient.cc
        EvCoefficient.cc
        FCoefficient.cc
        KappaCoefficient.cc
        UvCoefficient.cc)

set(libangcorr_src
        AngularCorrelation.cc
        W_dir_dir.cc
        W_pol_dir.cc)

set(libsampler_src
        AngCorrRejectionSampler.cc
        CascadeSampler.cc
        ReferenceFrameSampler.cc
        SphereIntegrator.cc
        SpherePointSampler.cc
        SphereRejectionSampler.cc
        SpotlightSampler.cc)

add_library(objcoeff OBJECT ${libcoeff_src})
add_library(objangcorr OBJECT ${libangcorr_src})
add_library(objsampler OBJECT ${libsampler_src})
target_link_libraries(objcoeff PUBLIC GSL::gsl GSL::gslcblas PRIVATE alpaca_options)
target_link_libraries(objangcorr PUBLIC objcoeff PRIVATE alpaca_options)
target_link_libraries(objsampler PUBLIC GSL::gsl GSL::gslcblas objcoeff objangcorr PRIVATE alpaca_options)

set_property(TARGET objcoeff PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET objangcorr PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET objsampler PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(coeff SHARED $<TARGET_OBJECTS:objcoeff>)
add_library(coeff_static STATIC $<TARGET_OBJECTS:objcoeff>)
set_target_properties(coeff_static PROPERTIES OUTPUT_NAME coeff)
target_link_libraries(coeff PUBLIC GSL::gsl GSL::gslcblas PRIVATE alpaca_options)
target_link_libraries(coeff_static PUBLIC GSL::gsl GSL::gslcblas PRIVATE alpaca_options)
target_include_directories(objcoeff PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
add_library(alpaca::coeff ALIAS coeff)

add_library(angcorr SHARED $<TARGET_OBJECTS:objangcorr>)
add_library(angcorr_static STATIC $<TARGET_OBJECTS:objangcorr>)
set_target_properties(angcorr_static PROPERTIES OUTPUT_NAME angcorr)
target_link_libraries(angcorr PUBLIC coeff PRIVATE alpaca_options)
target_link_libraries(angcorr_static PUBLIC coeff_static PRIVATE alpaca_options)
target_include_directories(objangcorr PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
add_library(alpaca::angcorr ALIAS angcorr)

add_library(sampler SHARED $<TARGET_OBJECTS:objsampler>)
add_library(sampler_static STATIC $<TARGET_OBJECTS:objsampler>)
set_target_properties(sampler_static PROPERTIES OUTPUT_NAME sampler)
target_link_libraries(sampler PUBLIC GSL::gsl GSL::gslcblas coeff angcorr PRIVATE alpaca_options)
target_link_libraries(sampler_static PUBLIC GSL::gsl GSL::gslcblas coeff_static angcorr_static PRIVATE alpaca_options)
target_include_directories(objsampler PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
add_library(alpaca::sampler ALIAS sampler)

if(BUILD_PYTHON)
        find_package(pybind11 QUIET)
        if(NOT pybind11_FOUND)
                include(FetchContent)
                FetchContent_Declare(
                        pybind11
                        GIT_REPOSITORY https://github.com/pybind/pybind11.git
                        GIT_TAG        v2.10.4)
                FetchContent_MakeAvailable(pybind11)
        endif()

        pybind11_add_module(alpaca PythonBindings.cc)
        target_link_libraries(alpaca PRIVATE coeff_static angcorr_static sampler_static alpaca_options)
        target_include_directories(
                alpaca PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                       $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
endif()
