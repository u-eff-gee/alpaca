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

if(BUILD_TESTS)
    include_directories(${PROJECT_SOURCE_DIR}/include/)

    enable_testing()

    add_executable(test_f_coefficient test_f_coefficient.cc)
    target_link_libraries(test_f_coefficient fCoefficient)
    add_test(test_f_coefficient test_f_coefficient)

    add_executable(test_uv_coefficient test_Uv_coefficient.cc)
    target_link_libraries(test_uv_coefficient uvCoefficient)
    add_test(test_uv_coefficient test_uv_coefficient)

    add_executable(test_w_dir_dir test_w_dir_dir.cc)
    target_link_libraries(test_w_dir_dir w_dir_dir)
    add_test(test_w_dir_dir test_w_dir_dir)

    add_executable(test_kappa_coefficient test_kappa_coefficient.cc)
    target_link_libraries(test_kappa_coefficient kappa_coefficient)
    add_test(test_kappa_coefficient test_kappa_coefficient)

    add_executable(test_ev_coefficient test_Ev_coefficient.cc)
    target_link_libraries(test_ev_coefficient alphavCoefficient evCoefficient)
    add_test(test_ev_coefficient test_ev_coefficient)

    add_executable(test_w_pol_dir test_w_pol_dir.cc)
    target_link_libraries(test_w_pol_dir w_pol_dir ${GSL_LIBRARIES})
    add_test(test_w_pol_dir test_w_pol_dir)

    add_executable(test_unobserved test_unobserved.cc)
    target_link_libraries(test_unobserved w_dir_dir ${GSL_LIBRARIES})
    add_test(test_unobserved test_unobserved)

    add_executable(test_euler_angle_rotation test_euler_angle_rotation.cc)
    target_link_libraries(test_euler_angle_rotation euler_angle_rotation w_dir_dir w_pol_dir ${GSL_LIBRARIES})
    add_test(test_euler_angle_rotation test_euler_angle_rotation)

    add_executable(test_angular_correlation test_angular_correlation.cc)
    target_link_libraries(test_angular_correlation angular_correlation euler_angle_rotation w_dir_dir w_pol_dir ${GSL_LIBRARIES})
    add_test(test_angular_correlation test_angular_correlation)

    add_executable(test_angular_correlation_io test_angular_correlation_io.cc)
    target_link_libraries(test_angular_correlation_io angular_correlation)
    add_test(test_angular_correlation_io test_angular_correlation_io)

    add_executable(test_elliptic_integral test_elliptic_integral.cc)
    target_link_libraries(test_elliptic_integral spherePointSampler ${GSL_LIBRARIES})
    add_test(test_elliptic_integral test_elliptic_integral)

    add_executable(test_spiral_segment test_spiral_segment.cc)
    target_link_libraries(test_spiral_segment spherePointSampler ${GSL_LIBRARIES})
    add_test(test_spiral_segment test_spiral_segment)

    add_executable(test_find_c test_find_c.cc)
    target_link_libraries(test_find_c spherePointSampler)
    add_test(test_find_c test_find_c)

    add_executable(test_sphere_integrator test_sphere_integrator.cc)
    target_link_libraries(test_sphere_integrator sphereIntegrator)
    add_test(test_sphere_integrator test_sphere_integrator)

    add_executable(test_normalization test_normalization.cc)
    target_link_libraries(test_normalization sphereIntegrator w_dir_dir w_pol_dir)
    add_test(test_normalization test_normalization)

    add_executable(test_find_Theta_j test_find_Theta_j.cc)
    target_link_libraries(test_find_Theta_j spherePointSampler)
    add_test(test_find_Theta_j test_find_Theta_j)

    add_executable(test_upper_limit test_upper_limit.cc)
    target_link_libraries(test_upper_limit spherePointSampler w_dir_dir w_gamma_gamma w_pol_dir)

    add_executable(test_sphere_rejection_sampler test_sphere_rejection_sampler.cc)
    target_link_libraries(test_sphere_rejection_sampler sphereRejectionSampler w_dir_dir)

    add_executable(test_string_representation test_string_representation.cc)
    target_link_libraries(test_string_representation w_dir_dir w_pol_dir)
    add_test(test_string_representation test_string_representation)

    add_executable(test_memory_leak test_memory_leak.cc)
    target_link_libraries(test_memory_leak angular_correlation)
    add_test(test_memory_leak test_memory_leak)
endif(BUILD_TESTS)