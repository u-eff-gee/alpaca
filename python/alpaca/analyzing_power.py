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

import warnings

import numpy as np

from alpaca.angular_correlation import AngularCorrelation
from alpaca.transition import Transition

CONVENTION = {"natural": 1.0, "KPZ": -1.0}


def arctan_grid(n, abs_delta_max=100.0):
    r"""Create an equidistant grid for the arctangent of the multipole-mixing ratio

    Multipole mixing ratios, as defined in the present code, are variables whose range includes
    the entire set of real numbers.
    This is because they are defined - as their name suggests - as a ratio of real-valued
    matrix elements.
    However, the largest variations of the analyzing power with the mixing ratio usually occur
    for values on the order of unity, which imply a strong competition between the two possible
    multipoles.
    Due to this expected behavior, it is not efficient to evaluate the analyzing power on an
    equidistant grid in the mixing ratio.
    Instead, this function samples :math:`N` equidistant values
    :math:`\mathrm{arctan} \left(\delta_i\right)`
    (:math:`i \in [0, N-1]`) in the range
    :math:`\left[ -\delta_\mathrm{max}, \delta_\mathrm{max}\right]`
    (the range includes both :math:`-\delta_\mathrm{max}` and :math:`-\delta_\mathrm{max}`)
    on an arctangent-compressed axis, and returns a list of the corresponding :math:`\delta_i`.
    The arctangent transformation does exactly what is needed, i.e. compressing the infinitely
    large ranges where the alternative multipole dominates, while keeping a high grid density
    in regions around 0 where the contributions by the two multipoles are approximately equal,
    or the primary multipole dominates.

    Parameters
    ----------
    n: int
        Number of grid points. Must be larger than 2 to at least include the two limits.
        Note that only odd numbers of points will include a mixing ratio of exactly zero.
    abs_delta_max: float
        Maximum absolute value of the multipole mixing ratio which determines the limits of
        the symmetric grid (default: 100).

    Returns
    -------
    ndarray
        :math:`\delta_i`, grid points for the multipole-mixing ratio.

    Warns
    -----
    UserWarning
        If the number of grid points is smaller than 2 or even.
    """

    if n < 2:
        warnings.warn(
            "The number of grid points must be larger than one, so that at least the two limits of the interval can be included. Using n=2."
        )
        n = 2
    if n % 2 == 0:
        warnings.warn(
            "An even number of grid points was given. While this is a perfectly valid input, please be aware that the set of points will not include a multipole mixing of exactly zero. Using any valid odd number will include it."
        )
    arctan_delta_max = np.arctan(np.abs(abs_delta_max))
    arctan_deltas = np.linspace(-arctan_delta_max, arctan_delta_max, n)
    return np.tan(arctan_deltas)


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

    def __call__(self, theta=0.5 * np.pi):
        r"""Evaluate the analyzing power

        Parameters
        ----------
        theta: float or ndarray
            Polar angle :math:`\theta` in radians (default: 90 degrees).

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

    def evaluate(self, delta, delta_values, theta=0.5 * np.pi):
        r"""Evaluate the analyzing power for a given value of the multipole mixing ratio

        Based on the cascade in this AnalyzingPower object, this function creates a new object
        with the given values of the multipole mixing ratio.
        It is assumed that only one variable is needed to obtain all the mixing ratios of the
        cascade.

        Parameters
        ----------
        delta: float or ndarray
            Value of the variable.
        delta_values: list of str or float or callable
            This list indicates which of the mixing ratios in the cascade are variables
            (arbitrary string), dependent on the variables (callable), or should be fixed to a
            value (float).
            The mixing-ratio values are given in the same order as the cascade transitions in
            `alpaca.AngularCorrelation`.
            Note that `alpaca.AnalyzingPower.evaluate` can only evaluate the analyzing power for
            a single variable, i.e. putting in different strings for different transitions will
            not make this a higher-dimensional function.
            The option of giving a callable function instead of a value can be used to simulate
            what would happen if the wrong mixing-ratio convention had been used.
            For example, in the Biedenharn convention, the mixing ratios for both mixing ratios
            in an elastic two-step cascade are the same: `delta_values=["delta", "delta"]`.
            In the Krane-Steffen-Wheeler convention, however, the first mixing ratio would have
            the opposite sign.
            To achieve this, put `delta_values=["delta", lambda x: -x]`.
        theta: float or ndarray
            Polar angle :math:`\theta` in radians (default: 90 degrees).

        Returns
        -------
        float
            Value of the analyzing power at the given multipole mixing ratio.
        """
        original_shape = np.shape(delta)
        scalar_output = isinstance(delta, (int, float))
        delta = np.reshape(delta, (np.size(delta),))
        asymmetries = np.zeros(len(delta))

        for i, d in enumerate(delta):
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
                                d,
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
                                delta_values[j](d),
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
            asymmetries[i] = AnalyzingPower(
                AngularCorrelation(
                    self.angular_correlation.initial_state, cascade_steps
                ),
                PQ=self.PQ,
                convention=self.convention,
            )(theta)
        if scalar_output:
            return asymmetries[0]
        return np.reshape(asymmetries, original_shape)
