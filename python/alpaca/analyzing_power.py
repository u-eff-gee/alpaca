# This file is part of alpaca.

# alpaca is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# alpaca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with alpaca.  If not, see <https://www.gnu.org/licenses/>.

# Copyright (C) 2021 Udo Friman-Gayer

import numpy as np

CONVENTION = {"natural": 1.0, "KPZ": -1.0}


class AnalyzingPower:
    r"""Analyzing power for a given angular correlation

    In general, the analyzing power :math:`A` relates the angular correlation
    :math:`W \left( \theta, \varphi \right)` in two different directions
    :math:`\left( \theta, \varphi \right)` and
    :math:`\left( \theta^\prime, \varphi^\prime \right)` as

    ..math:: \frac{W \left( \theta, \varphi \right) - W \left( \theta^\prime, \varphi^\prime \right)}{W \left( \theta, \varphi \right) + W \left( \theta^\prime, \varphi^\prime \right)},

    i.e. it quantifies the relative difference of the probabilities.

    Due to the symmetry of the angular correlations in nuclear physics experiments, they often use
    a definition where :math:`\theta = \theta^\prime` and :math:`\varphi, \varphi^\prime \in \left\{ 0, \pi/2 \right\}`.
    The most recent review article on the nuclear resonance fluorescence technique by Kneissl,
    Pitz, and Zilges (KPZ) :cite:`Kneissl1996` gives the analyzing power as:

    ..math:: A \left( \theta \right) = \frac{W \left( \theta, \varphi = \pi/2 \right) - W \left( \theta, \varphi = 0 \right)}{W \left( \theta, \varphi = \pi/2 \right) + W \left( \theta, \varphi = 0 \right)}

    One frequently also finds another definition of the analyzing power that differs from the KPZ
    convention by a minus sign:

    ..math:: A \left( \theta \right) = \frac{W \left( \theta, \varphi = 0 \right) - W \left( \theta, \varphi = \pi / 2 \right)}{W \left( \theta, \varphi = \pi/2 \right) + W \left( \theta, \varphi = 0 \right)}

    The latter convention may originate from the investigation of the dipole response of even-even
    nuclei with polarized photon beams, as pioneered by Pietralla et al. :cite:`Pietralla2001`.
    First experiments of this kind used detectors in a cross-shaped arrangement
    (:math:`\theta = \pi/2`).
    In this case, the analyzing powers for :math:`0^+ \to 1^+ \to 0^+` and
    :math:`0^+ \to 1^- \to 0^+` conveniently evaluate to :math:`A \left( \pi/2 \right) = +1`
    and :math:`A \left( \pi/2 \right) = -1`, respectively.

    This class supports both the 'KPZ' and the 'natural' convention.

    Attributes
    ----------
    angular_correlation: AngularCorrelation object
        Angular correlation
    convention: str
        'KPZ' for the convention of Kneissl, Pitz and Zilges, or 'natural' (default) for the natural convention.
    """

    def __init__(self, angular_correlation, convention="natural"):
        """Initialization

        Define the angular correlation for which the analyzing power should be evaluated and the
        convention.

        Parameters
        ----------
        angular_correlation: AngularCorrelation object
            Angular correlation
        convention: str
            'KPZ' for the convention of Kneissl, Pitz and Zilges, or 'natural' (default) for the natural convention.
        """
        self.angular_correlation = angular_correlation
        self.convention = convention

    def __call__(self, theta):
        r"""Evaluate the analyzing power

        Parameters
        ----------
        theta: float or ndarray
            Polar angle :math:`\theta` in radians.

        Returns
        -------
        float, :math:`A \left( \theta \right)`
        """
        w_para = self.angular_correlation(theta, 0.0)
        w_perp = self.angular_correlation(theta, 0.5 * np.pi)

        return CONVENTION[self.convention] * (w_para - w_perp) / (w_para + w_perp)
