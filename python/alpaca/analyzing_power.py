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
    r"""Find the intersection of two intervals

    Given two intervals :math:`\left[ a, b \right]` and :math:`\left[ c, d \right]`, this
    function finds their intersection.

    When searching for the intersection, all comparisons are made with the :math:`\leq` and
    :math:`\geq` operators.
    This means that if both arrays only have a single real number in common, i.e.
    :math:`\left[ a, b \right]` and :math:`\left[ b, d \right]`, the returned intersection will
    be :math:`\left[ b, b \right]`, and not an empty list.

    The intervals are sorted before they are processed by this function, so the interval limits can
    be unsorted.

    Parameters
    ----------
    interval_1: list or ndarray of two float
        First interval. May be unsorted.
    interval_2: list or ndarray of two float
        Second interval. May be unsorted.

    Returns
    -------
    list of two float or empty list
        Intersection of the two intervals as a new interval, or an empty list.

    Examples
    --------
    >>> intersection_of_two_intervals([0.0, 1.0], [0.5, 1.5])
    [0.5, 1.0]
    >>> intersection_of_two_intervals([0.0, 1.0], [0.5, 0.6])
    [0.5, 0.6]
    >>> intersection_of_two_intervals([0.0, 1.0], [1.0, 1.5])
    [1.0, 1.0]
    >>> intersection_of_two_intervals([0.0, 1.0], [1.5, 2.0])
    []
    """
    intersection = []

    interval_1 = np.sort(interval_1)
    interval_2 = np.sort(interval_2)

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
    r"""Find the intersection of an interval with a list of intervals

    Given an interval :math:`\left[ a, b \right]` and a set of intervals :math:`\left\{ \left[ c, d \right], \left[ e, f \right], ... \right \}`, this functions determines the intersections of the former with all elements of the latter.

    See also `alpaca.analyzing_power.intersection_of_two_intervals`.

    Parameters
    ----------
    interval_1: list or ndarray of two float
        First interval. May be unsorted.
    list_of_intervals: list of lists or ndarrays of two float
        List of intervals. The single intervals may be unsorted.

    Returns
    -------
    list of lists of two float, or empty list
        Intersection of the interval with the list of intervals as a new list of intervals, or an empty list.

    Examples
    --------
    >>> intersection_of_interval_with_list_of_intervals([0.0, 1.0], [[0.0, 0.5], [0.8, 1.5]])
    [[0.0, 0.5], [0.8, 1.0]]
    """
    intersections = []

    for interval in list_of_intervals:
        intersection = intersection_of_two_intervals(interval_1, interval)
        if len(intersection) > 0:
            intersections.append(intersection)

    return intersections


def intersection(list_of_intervals_1, list_of_intervals_2):
    r"""Find the intersections of a list of intervals with another list of intervals

    Given a set of intervals :math:`\left\{ \left[ \alpha, \beta \right], \left[ \gamma, \delta \right], ... \right \}` and a set of intervals :math:`\left\{ \left[ a, b \right], \left[ c, d \right], ... \right \}`, this functions determines the intersections of all of the former intervals with all elements of the latter.

    See also `alpaca.analyzing_power.intersection_of_two_intervals` and `alpaca.analyzing_power.intersection_of_interval_with_list_of_intervals`.

    Parameters
    ----------
    list_of_intervals_1: list of lists or ndarrays of two float
        First list of intervals. The single intervals may be unsorted.
    list_of_intervals_2: list of lists or ndarrays of two float
        Second list of intervals. The single intervals may be unsorted.

    Returns
    -------
    list of lists of two float, or empty list
        Intersection of the first list of intervals with the second list of intervals as a new list of intervals, or an empty list.

    Examples
    --------
    >>> intersection([[0.0, 0.1], [0.3, 0.4], [0.6, 1.0]], [[0.3, 0.5], [0.8, 0.9]])
    [[0.3, 0.4], [0.8, 0.9]]
    """
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
        asymmetry,
        delta_values,
        theta=0.5 * np.pi,
        n_delta=int(1e3),
        atol=1e-3,
        abs_delta_max=100.0,
        return_intervals=False,
    ):
        r"""Find ranges of multipole mixing ratios that are consistent with a given asymmetry

        A common task in nuclear spectroscopy is to determine the multipole mixing ratio of a
        transition from an observed angular correlation, which is manifested in a count-rate
        asymmetry between two detectors :math:`\epsilon_\mathrm{exp}`.
        The precision of the mixing-ratio measurement does not only depend on the precision of
        the asymmetry, given as a coverage interval :math:`\left[ \epsilon_\mathrm{exp}^\mathrm{min}, \epsilon_\mathrm{exp}^\mathrm{max}\right]`, but also on the cascade of interest.
        In general, multiple values of the mixing ratio can lead to the same asymmetry.

        This function assumes that a cascade is studied that has a single unknown mixing ratio
        :math:`delta`.
        In general, :math:`delta` may appear in more than one transition of the cascade.
        The simplest example is an 'elastic' process with a ground state - excited state - ground
        state sequence.
        In order to solve the system of inequalities

        ..math:: \epsilon_\mathrm{exp}^\mathrm{min} \leq P Q A \left( \delta \right) \leq \epsilon_\mathrm{exp}^\mathrm{max},

        this function evaluates :math:`\epsilon`/:math:`A` on an equidistant grid in
        `\mathrm{arctan} \left( \delta_i \right)`.
        The index :math:`i` is restricted to :math:`0 < i < N-1`, where :math:`N` is
        the number of grid points.
        In the given set of :math:`\epsilon_i = P Q A \left( \delta_ i \right)`, the function then finds the values of :math:`\delta_i` that fulfil the
        inequalities.
        If several `\delta_i` in a row fulfil the inequalities, they may belong to a contiguous
        interval.
        This function either returns a list of `\delta_i`, or these presumedly contiguous intervals.
        Obviously, a high enough value for :math:`N` must be chosen to keep numerical uncertainties
        at a minimum compared to the experimental value.

        In an alternative mode, the asymmetry can also be given as a single value :math:`\epsilon_\mathrm{exp}`.
        In this case, the inequality

        ..math:: \left| \epsilon_\mathrm{exp} - P Q A \left( \delta \right) \right| \leq t,

        is evaluated, where :math:`t` is a user-defined tolerance.
        The second mode requires :math:`N` and :math:`t` to be well matched.

        Parameters
        ----------
        asymmetry: float or [float, float]
            Value of the asymmetry or interval of possible values (coverage interval,
            standard deviations, ...)
        delta_values: list of str or float or callable
            This list indicates which of the mixing ratios in the cascade are variables
            (arbitrary string), dependent on the variables (callable), or should be fixed to a 
            value (float).
            The mixing-ratio values are given in the same order as the cascade transitions in `alpaca.AngularCorrelation`.
            Note that `alpaca.AnalyzingPower.find_delta_brute_force` can only solve for a single
            variable, i.e. putting in different strings for different transitions will not make
            this function solve a multi-parameter problem.
            The option of giving a callable function instead of a value can be used to simulate 
            what would happen if the wrong mixing-ratio convention had been used.
            For example, in the Biedenharn convention, the mixing ratios for both mixing ratios
            in an elastic two-step cascade are the same: `delta_values=["delta", "delta"]`.
            In the Krane-Steffen-Wheeler convention, however, the second mixing ratio would have
            a different sign.
            To achieve this, put `delta_values=["delta", lambda x: -x]`.
        theta: float or ndarray
            Polar angle :math:`\theta` in radians (default: 90 degrees).
        n_delta: int
            :math:`N`, number of trial multipole-mixing ratios (default: 1000).
        atol: float
            :math:`t`, absolute tolerance for determining the numerical equality of a calculate
            value and the given single-value asymmetry (default: 0.001).
        abs_delta_max: float
            Maximum value of the mixing ratio, which determines the positive and negative limits
            of the search grid (default: 100).
        return_intervals: bool
            Determines whether a list of :math:`\delta_i` (False) or a list of intervals (True)
            should be returned (default: False).

        Returns
        -------
        list of float or list of [float, float]
            Depending on the `return_intervals` setting, returns a list of all :math:`\delta_i`
            that match the given asymmetry (interval) or a list of intervals of matching values.
        """
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
                elif callable(delta_values[j]):
                    cascade_steps.append(
                        [
                            Transition(
                                cas_ste[0].em_char,
                                cas_ste[0].two_L,
                                cas_ste[0].em_charp,
                                cas_ste[0].two_Lp,
                                delta_values[j](delta),
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

            if isinstance(asymmetry, (int, float)):
                if np.abs(ana_pow - asymmetry) < atol:
                    delta_results.append(delta)
                    delta_matches[i] = True
            else:
                if asymmetry[0] - atol <= ana_pow <= asymmetry[1] + atol:
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
