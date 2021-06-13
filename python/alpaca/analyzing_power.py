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

from alpaca.angular_correlation import AngularCorrelation
from alpaca.transition import Transition

CONVENTION = {"natural": 1.0, "KPZ": -1.0}


def intersection_of_two_intervals(interval_1, interval_2):
    intersection = []

    if interval_1[0] > interval_2[1]:
        return []
    elif interval_1[0] >= interval_2[0] and interval_1[0] <= interval_2[1]:
        if interval_1[1] <= interval_2[1]:
            return [interval_1[0], interval_1[1]]
        else:
            return [interval_1[0], interval_2[1]]
    elif interval_1[1] >= interval_2[0]:
        if interval_1[1] <= interval_2[1]:
            return [interval_2[0], interval_1[1]]
        else:
            return [interval_2[0], interval_2[1]]

    return []


def intersection_of_interval_with_list_of_intervals(interval_1, list_of_intervals):
    intersections = []

    for interval in list_of_intervals:
        intersection = intersection_of_two_intervals(interval_1, interval)
        if len(intersection) > 0:
            intersections.append(intersection)

    return intersections


def intersection(list_of_intervals_1, list_of_intervals_2):
    intersections = []

    for interval_1 in list_of_intervals_1:
        intersection = intersection_of_interval_with_list_of_intervals(
            interval_1, list_of_intervals_2
        )
        if len(intersection) > 0:
            for interval in intersection:
                intersections.append(interval)

    return intersections


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

    The analyzing power quantifies the count-rate ratio between two detectors that form a
    'polarimeter' in the limit where all geometric effects (target size, detector solid angle,
    attenuation, ...) can be neglected, and the incoming photon beam has a perfect linear
    polarization.
    In a real experiment, the beam polarization :math:`P` and the so-called polarization
    sensitivity :math:`Q` of the polarimeter will lead to an observed asymmetry :math:`\epsilon`
    that is related to the analyzing power via :cite:`Kneissl1996`:

    ..math:: \epsilon \left( \theta \right) = P Q A \left( \theta \right).

    By rotation of the polarization axis or a superposition multiple linearly polarized photon
    beams, the polarization can vary in the range

    ..math:: P \in \left[ -1, 1\right].

    The default value of the polarization is :math:`P = 1`, corresponding to an electric-field
    vector along the :math:`x` axis only.
    A value of :math:`P = -1` corresponds to a rotation of the :math:`P = 1` polarization plane by
    90 degrees.

    The polarization sensitivity :math:`Q` varies within the range

    ..math:: Q \in \left[ 0, 1 \right],

    i.e. it decreases the analyzing power.
    It should be noted that :math:`Q` is in general different for each angular correlation.

    The product :math:`PQ` can be given to this function as an additional parameter to return an
    asymmetry instead of an analyzing power.

    Attributes
    ----------
    angular_correlation: AngularCorrelation object
        Angular correlation
    PQ: float
        Product of the photon polarization and the polarization sensitivity.
    convention: str
        'KPZ' for the convention of Kneissl, Pitz and Zilges, or 'natural' (default) for the natural convention.

    Raises
    ------
    ValueError
        If the absolute value of PQ is larger than 1.
    """

    def __init__(self, angular_correlation, PQ=1.0, convention="natural"):
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
        if np.abs(PQ) > 1.0:
            raise ValueError(
                "The absolute value of PQ must be smaller than or equal to 1."
            )
        self.PQ = PQ
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

        return (
            self.PQ
            * CONVENTION[self.convention]
            * (w_para - w_perp)
            / (w_para + w_perp)
        )

    def find_delta_brute_force(
        self,
        analyzing_power_value,
        delta_values,
        theta,
        n_delta=int(1e3),
        atol=1e-3,
        abs_delta_max=100.0,
        return_intervals=False,
    ):

        arctan_delta_max = np.arctan(abs_delta_max)
        arctan_deltas = np.linspace(-arctan_delta_max, arctan_delta_max, n_delta)
        deltas = np.tan(arctan_deltas)

        delta_results = []
        delta_matches = [False] * n_delta

        for i, delta in enumerate(deltas):
            cascade_steps = []
            for j, cas_ste in enumerate(self.angular_correlation.cascade_steps):
                if isinstance(delta_values[j], str):
                    cascade_steps.append(
                        [
                            Transition(
                                cas_ste[0].em_char,
                                cas_ste[0].two_L,
                                cas_ste[0].em_charp,
                                cas_ste[0].two_Lp,
                                delta,
                            ),
                            cas_ste[1],
                        ]
                    )
                else:
                    cascade_steps.append(
                        [
                            Transition(
                                cas_ste[0].em_char,
                                cas_ste[0].two_L,
                                cas_ste[0].em_charp,
                                cas_ste[0].two_Lp,
                                delta_values[j],
                            ),
                            cas_ste[1],
                        ]
                    )
            ana_pow = AnalyzingPower(
                AngularCorrelation(
                    self.angular_correlation.initial_state, cascade_steps
                ),
                convention=self.convention,
            )(theta)

            if isinstance(analyzing_power_value, (int, float)):
                if np.abs(ana_pow - analyzing_power_value) < atol:
                    delta_results.append(delta)
                    delta_matches[i] = True
            else:
                if (
                    analyzing_power_value[0] - atol
                    <= ana_pow
                    <= analyzing_power_value[1] + atol
                ):
                    delta_results.append(delta)
                    delta_matches[i] = True

        if return_intervals:
            intervals = []
            interval_start = None
            for i, matching in enumerate(delta_matches):
                if matching and interval_start is None:
                    interval_start = deltas[i]
                if not matching and interval_start is not None:
                    intervals.append([interval_start, deltas[i]])
                    interval_start = None
                if i == n_delta - 1 and interval_start is not None:
                    intervals.append([interval_start, deltas[-1]])
            return intervals

        return (delta_results, delta_matches)
