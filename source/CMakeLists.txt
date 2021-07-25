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

include_directories(${PROJECT_SOURCE_DIR}/include/)

add_library(fCoefficient FCoefficient.cc)
target_link_libraries(fCoefficient ${GSL_LIBRARIES})

add_library(avCoefficient AvCoefficient.cc)
target_link_libraries(avCoefficient fCoefficient)

add_library(w_gamma_gamma W_gamma_gamma.cc)
target_compile_options(w_gamma_gamma PUBLIC -Wall -Wextra)

add_library(w_dir_dir W_dir_dir.cc)
target_link_libraries(w_dir_dir avCoefficient uvCoefficient)

add_library(kappa_coefficient KappaCoefficient.cc)
target_link_libraries(kappa_coefficient ${GSL_LIBRARIES})

add_library(alphavCoefficient AlphavCoefficient.cc)
target_link_libraries(alphavCoefficient fCoefficient kappa_coefficient)

add_library(evCoefficient EvCoefficient.cc)
target_link_libraries(evCoefficient fCoefficient)

add_library(uvCoefficient UvCoefficient.cc)
target_link_libraries(uvCoefficient ${GSL_LIBRARIES})

add_library(w_pol_dir W_pol_dir.cc)
target_link_libraries(w_pol_dir alphavCoefficient avCoefficient w_dir_dir)

add_library(angular_correlation SHARED AngularCorrelation.cc)
target_link_libraries(angular_correlation euler_angle_rotation w_dir_dir w_pol_dir)

add_library(euler_angle_rotation EulerAngleRotation.cc)

add_library(spherePointSampler SpherePointSampler.cc)
target_link_libraries(spherePointSampler ${GSL_LIBRARIES})

add_library(sphereIntegrator SphereIntegrator.cc)
target_link_libraries(sphereIntegrator spherePointSampler)

add_library(sphereRejectionSampler SphereRejectionSampler.cc)

add_library(angcorrRejectionSampler AngCorrRejectionSampler.cc)
target_link_libraries(angcorrRejectionSampler sphereRejectionSampler angular_correlation)