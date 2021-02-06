configure_file(setup.py setup.py)
configure_file(setup.cfg setup.cfg)
configure_file(tox.ini tox.ini)

configure_file(alpaca/analyzing_power.py alpaca/analyzing_power.py)
configure_file(alpaca/angular_correlation.py alpaca/angular_correlation.py @ONLY)
configure_file(alpaca/level_scheme_plotter.py alpaca/level_scheme_plotter.py)
configure_file(alpaca/__init__.py alpaca/__init__.py)
configure_file(alpaca/state.py alpaca/state.py)
configure_file(alpaca/transition.py alpaca/transition.py)

configure_file(test/__init__.py test/__init__.py)
configure_file(test/test_io.py test/test_io.py)
configure_file(test/test_level_scheme_plotter.py test/test_level_scheme_plotter.py)