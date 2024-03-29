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

configure_file(setup.py setup.py)
configure_file(setup.cfg setup.cfg)
configure_file(tox.ini.in tox.ini)

configure_file(alpaca/analyzing_power.py alpaca/analyzing_power.py)
configure_file(alpaca/analyzing_power_plotter.py alpaca/analyzing_power_plotter.py)
configure_file(alpaca/angular_correlation.py alpaca/angular_correlation.py @ONLY)
configure_file(alpaca/angular_correlation_plotter.py alpaca/angular_correlation_plotter.py)
configure_file(alpaca/angular_correlation_table.py alpaca/angular_correlation_table.py)
configure_file(alpaca/interval_intersections.py alpaca/interval_intersections.py)
configure_file(alpaca/inversion_by_grid_evaluation.py alpaca/inversion_by_grid_evaluation.py)
configure_file(alpaca/inversion_by_piecewise_interpolation.py alpaca/inversion_by_piecewise_interpolation.py)
configure_file(alpaca/level_scheme_plotter.py alpaca/level_scheme_plotter.py)
configure_file(alpaca/__init__.py alpaca/__init__.py)
configure_file(alpaca/state.py alpaca/state.py)
configure_file(alpaca/transition.py alpaca/transition.py)

configure_file(test/__init__.py test/__init__.py)
configure_file(test/test_analyzing_power.py test/test_analyzing_power.py)
configure_file(test/test_analyzing_power_pure_transitions.py test/test_analyzing_power_pure_transitions.py)
configure_file(test/test_analyzing_power_single_mixing.py test/test_analyzing_power_single_mixing.py)
configure_file(test/test_angular_correlation.py test/test_angular_correlation.py)
configure_file(test/test_angular_correlation_plotter.py test/test_angular_correlation_plotter.py)
configure_file(test/test_angular_correlation_table.py test/test_angular_correlation_table.py)
configure_file(test/test_interval_intersections.py test/test_interval_intersections.py)
configure_file(test/test_inversion_by_grid_evaluation.py test/test_inversion_by_grid_evaluation.py)
configure_file(test/test_inversion_by_piecewise_interpolation.py test/test_inversion_by_piecewise_interpolation.py)
configure_file(test/test_io.py test/test_io.py)
configure_file(test/test_level_scheme_plotter.py test/test_level_scheme_plotter.py)
configure_file(test/test_memory_leak.py test/test_memory_leak.py)
configure_file(test/test_state.py test/test_state.py)
configure_file(test/test_transition.py test/test_transition.py)