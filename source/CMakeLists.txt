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

add_library(fCoefficient FCoefficient.cc)
target_link_libraries(fCoefficient ${GSL_LIBRARIES})
target_include_directories(fCoefficient PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(fCoefficient PROPERTIES PUBLIC_HEADER include/alpaca/FCoefficient.hh)

add_library(avCoefficient AvCoefficient.cc)
target_link_libraries(avCoefficient fCoefficient)
target_include_directories(avCoefficient PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(avCoefficient PROPERTIES PUBLIC_HEADER include/alpaca/AvCoefficient.hh)

add_library(w_dir_dir W_dir_dir.cc)
target_link_libraries(w_dir_dir avCoefficient uvCoefficient)
target_include_directories(w_dir_dir PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(w_dir_dir PROPERTIES PUBLIC_HEADER include/alpaca/W_dir_dir.hh)

add_library(kappa_coefficient KappaCoefficient.cc)
target_link_libraries(kappa_coefficient ${GSL_LIBRARIES})
target_include_directories(kappa_coefficient PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(kappa_coefficient PROPERTIES PUBLIC_HEADER include/alpaca/KappaCoefficient.hh)

add_library(alphavCoefficient AlphavCoefficient.cc)
target_link_libraries(alphavCoefficient fCoefficient kappa_coefficient)
target_include_directories(alphavCoefficient PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(alphavCoefficient PROPERTIES PUBLIC_HEADER include/alpaca/AlphavCoefficient.hh)

add_library(evCoefficient EvCoefficient.cc)
target_link_libraries(evCoefficient fCoefficient)
target_include_directories(evCoefficient PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(evCoefficient PROPERTIES PUBLIC_HEADER include/alpaca/EvCoefficient.hh)

add_library(uvCoefficient UvCoefficient.cc)
target_link_libraries(uvCoefficient ${GSL_LIBRARIES})
target_include_directories(uvCoefficient PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(uvCoefficient PROPERTIES PUBLIC_HEADER include/alpaca/UvCoefficient.hh)

add_library(w_pol_dir W_pol_dir.cc)
target_link_libraries(w_pol_dir alphavCoefficient avCoefficient w_dir_dir)
target_include_directories(w_pol_dir PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(w_pol_dir PROPERTIES PUBLIC_HEADER include/alpaca/W_pol_dir.hh)

add_library(angular_correlation SHARED AngularCorrelation.cc)
target_link_libraries(angular_correlation PUBLIC w_dir_dir w_pol_dir)
target_include_directories(angular_correlation PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(angular_correlation PROPERTIES PUBLIC_HEADER include/alpaca/AngularCorrelation.hh)

add_library(angular_correlation_static STATIC AngularCorrelation.cc)
target_link_libraries(angular_correlation_static PUBLIC w_dir_dir w_pol_dir)
target_include_directories(
  angular_correlation_static
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(
  angular_correlation_static PROPERTIES PUBLIC_HEADER
                                 include/alpaca/AngularCorrelation.hh)
set_target_properties(angular_correlation_static PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(spherePointSampler SpherePointSampler.cc)
target_link_libraries(spherePointSampler ${GSL_LIBRARIES})
target_include_directories(spherePointSampler PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(spherePointSampler PROPERTIES PUBLIC_HEADER include/alpaca/SpherePointSampler.hh)

add_library(sphereIntegrator SphereIntegrator.cc)
target_link_libraries(sphereIntegrator spherePointSampler)
target_include_directories(sphereIntegrator PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(sphereIntegrator PROPERTIES PUBLIC_HEADER include/alpaca/SphereIntegrator.hh)

add_library(referenceFrameSampler ReferenceFrameSampler.cc)
target_link_libraries(referenceFrameSampler)
target_include_directories(referenceFrameSampler PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(referenceFrameSampler PROPERTIES PUBLIC_HEADER include/alpaca/ReferenceFrameSampler.hh)

add_library(sphereRejectionSampler SphereRejectionSampler.cc)
target_link_libraries(sphereRejectionSampler referenceFrameSampler)
target_include_directories(sphereRejectionSampler PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(sphereRejectionSampler PROPERTIES PUBLIC_HEADER include/alpaca/SphereRejectionSampler.hh)

add_library(spotlightSampler SpotlightSampler.cc)
target_link_libraries(spotlightSampler referenceFrameSampler)
target_include_directories(sphereRejectionSampler PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_library(angcorrRejectionSampler AngCorrRejectionSampler.cc)
target_link_libraries(angcorrRejectionSampler sphereRejectionSampler angular_correlation)
target_include_directories(angcorrRejectionSampler PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(angcorrRejectionSampler PROPERTIES PUBLIC_HEADER include/alpaca/AngCorrRejectionSampler.hh)

add_library(cascadeSampler CascadeSampler.cc)
target_link_libraries(cascadeSampler angular_correlation angcorrRejectionSampler)
target_include_directories(cascadeSampler PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(cascadeSampler PROPERTIES PUBLIC_HEADER include/alpaca/CascadeSampler.hh)


if(BUILD_PYTHON)
        nanobind_add_module(alpaca alpaca.cc)
        target_link_libraries(alpaca PRIVATE angular_correlation_static)
        target_include_directories(
          alpaca PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
endif()
