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

from setuptools import setup, find_packages

setup(
        name='alpaca',
        version='@PROJECT_VERSION_MAJOR@.@PROJECT_VERSION_MINOR@.@PROJECT_VERSION_PATCH@',
        description='alpaca: a linearly-polarized angular-correlation application',
        url='http://github.com/uga-uga/alpaca',
        author='Udo Friman-Gayer',
        author_email='ufg@email.unc.edu',
        license='GPLv3',
        python_requires='>=3',
        packages=find_packages(),
        install_requires=['matplotlib', 'numpy'],
        setup_requires=['pytest-runner'],
        tests_require=['pytest', 'pytest-cov', 'tox']
)