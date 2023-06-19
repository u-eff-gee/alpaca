/*
    This file is part of alpaca.

    alpaca is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    alpaca is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with alpaca.  If not, see <https://www.gnu.org/licenses/>.

    Copyright (C) 2021-2023 Udo Friman-Gayer
*/

#pragma once

namespace alpaca {

/**
 * \brief Class for sampling the initial four-momenta of photons emitted in a
 * gamma-ray cascade.
 *
 * # Introduction
 *
 * The angular-correlation mechanism implemented in the alpaca code does not
 * only have an impact on the direction of emission of gamma rays in a cascade,
 * but also on their energy. Even if the level widths are usually on the order
 * of less than a few eV \cite Zilges2022, two related effects which are
 * introduced below can result in energy shifts on the order of a few keV. These
 * shifts are substantial enough to be resolved by state-of-the-art
 * high-resolution detectors.
 * Historically, these energy shifts delayed the observation of nuclear
 * resonance fluorescence, because they prevent the re-capture of gamma rays
 * from a radioactive source by another nucleus of the same isotope. The problem
 * can be overcome by the Mössbauer effect \cite Moessbauer1958 or by using
 * artificial photon sources.
 *
 * In the following, the two effects that cause an energy shift will be
 * introduced. Firstly, photoabsorption and -emission processes of photons in
 * the MeV energy range transfer a significant recoil energy to the
 * absorbing/emitting nucleus. The correction to the energy of the
 * absorbed/emitted photon is on the order of \cite Zilges2022
 *
 * \f[
 *      \frac{\Delta E}{2 M c^2} \Delta E < \approx 10^{-3} \times \Delta E,
 * \f]
 *
 * where \f$\Delta E\f$ is the energy difference between nuclear levels, and
 * \f$M\f$ is the mass of the nucleus. In order to obtain the numerical value of
 * \f$10^{-3}\f$, the fact that nuclear excited states have energies on the MeV
 * scale, while nuclear masses are on the GeV scale, has been used.
 *
 * Secondly, radiation absorbed or emitted by a nucleus in motion will be
 * subject to the Doppler effect. For a nucleus that recoils after absorbing a
 * photon, the energy shift is again on the order of (see the documentation of
 * the corresponding class methods)
 *
 * \f[
 *      \left( \frac{\sqrt{\left(M c^2 + \Delta E \right)^2 + E_\gamma^2}}{M
 * c^2} - 1 \right) E_\gamma < \approx 10^{-3} \times E_\gamma. \f]
 *
 * Here, \f$E_\gamma\f$ denotes the energy of the photon that was captured to
 * set the nucleus in motion. As seen above, \f$E_\gamma\f$ is approximately
 * equal to the level-energy difference up to the aforementioned recoil
 * correction. The same energy considerations as in the previous numerical
 * estimate were used here.
 *
 * # Theory
 *
 * Obviously, the distribution of possible momenta of gamma rays emitted in the
 * final transition of a cascade depends on the directions of emission of all
 * previous gamma rays. The FourMomentumSampler class uses a Monte Carlo
 * procedure to sample the four momenta of the emitted gamma rays, taking the
 * aforementioned recoil and Doppler effects into account.
 *
 * ## Assumptions
 *
 * The following assumptions are made for the sampling:
 *
 * - The lifetime of the states in the cascade is much shorter than the stopping
 * time of the moving atom in the environment. In other words, the nucleus moves
 * like a free particle. This assumption is often violated in reality. To
 * illustrate this, consider the realistic example of a popular NRF reaction:
 * After capturing a 4.4-MeV photon (the energy of its first excited state), a
 * \f$^{12}C\f$ nucleus moves at a velocity of about \f$0.03~c \approx 9 \times
 * 10^{-6}~\mathrm{mm}~\mathrm{fs}^{-1}\f$, which corresponds to a kinetic
 * energy of about \f$E_\gamma / 2 \approx 2.2 MeV\f$. A free nucleus would
 * deexcite by emitting a gamma ray in this state of motion. However, in
 * realistic conditions, the excited ion will interact with the atomic
 * environment and gradually lose its kinetic energy in collisions. If the
 * stopping time is comparable to or faster than the lifetime of the excited
 * state, the emission will occur at a lower velocity. In the present example,
 * the stopping power of a 2.2-MeV carbon ion in graphite (density: \f$2.2530 g
 * cm^{-3}\f$) is about \f$1.6 MeV mm^{-1}\f$ \cite Ziegler2007 \cite SRIM2022.
 * During the lifetime of the 4.4-MeV state of about 61 fs, the ion travels
 * about 0.0006 mm, corresponding to an energy loss of on the order of 1 keV. In
 * this case, the stopping would be negligible for practical purposes. However,
 * many lifetimes of nuclear states are orders of magnitude higher, and the
 * stopping becomes more effective for ions with higher proton numbers.
 * Therefore, the velocity of the nucleus may change significantly before the
 * deexcitation occurs.
 *
 * - The nucleus is at rest at the beginning of the reaction.
 * At usual laboratory temperatures of a few hundred Kelvin, this is well
 * justified, because the kinetic energy of even the lightest nuclei is only on
 * the order of a few dozen eV.
 *
 * ## Nomenclature
 *
 * Starting from these two assumptions, the absorption of a photon by a nucleus
 * and the subsequent emission of \f$n-1\f$ gamma rays can be described as an
 * \f$n\f$-step process. The following symbols will be used, some of which have
 * already been introduced.
 *
 * - \f$E_i\f$: Resonance energy of the \f$i\f$-th excited state of the nucleus.
 * The ground state is denoted as \f$i=0\f$.
 * - \f$E_\gamma\f$, \f$\mathbf{p}_\gamma\f$, \f$p_\gamma\f$: Energy and
 * momentum, and four-momentum of a photon.
 * - \f$E\f$, \f$\mathbf{p}\f$, \f$p\f$, \f$\beta\f$, \f$\gamma\f$: Energy,
 * momentum, four-momentum, and relativistic beta- and gamma-factor of the
 * nucleus.
 * - \f$M\f$: Mass of atom/ion/nucleus. The difference between the three masses
 * is small enough to be negligible for the present application. In fact, the
 * author is not sure which of the different possibilities would have to be used
 * in which step of the process. It is probably complicated, because the charge
 * state of the ion can change as it propagates through matter.
 * - \f$\theta\f$: Scattering angle of the photon in the laboratory frame with
 * respect to the previous step, which might be either an absorption or emission
 * process.
 *
 * The number of the step in the cascade will be indicated by prime symbols.
 * An index \f$R\f$ indicates that a quantity is measured in the rest frame of
 * the nucleus. All other quantities are assumed to be measured in the
 * laboratory frame.
 *
 * ## Capture
 *
 * In accordance with the assumptions above, start from a situation in which the
 * nucleus is at rest, and a photon is incoming. Without loss of generality, the
 * photon is assumed to propagate along the positive \f$z\f$ axis.
 *
 * \f[
 *      p_\gamma c = \left( E_\gamma, 0, 0, E_\gamma \right)
 * \f]
 * \f[
 *      p c = \left( M c^2, 0, 0, 0 \right)
 * \f]
 *
 * After the absorption of the photon, the compound nucleus is in an excited
 * state \f$i\f$, and traveling along the positive \f$z\f$ axis due to the
 * conservation of momentum ("recoil"):
 *
 * \f[
 *      p^\prime c = \left( \sqrt{\left[ M c^2 + \left( E_i - E_0 \right)
 * \right]^2 + |\mathbf{p}^\prime|^2 c^2}, \mathbf{p}^\prime c \right) \f]
 *
 * Conservation of four-momentum in the absorption process
 *
 * \f[
 *      p_\gamma + p = p^\prime
 * \f]
 *
 *  requires that
 *
 * \f[
 *      \mathbf{p}^\prime c = \left( 0, 0, E_\gamma \right)
 * \f]
 *
 * and -using the conservation of momentum in the expression for the
 * conservation of energy-
 *
 * \f[
 *      E_\gamma + M c^2 = \sqrt{\left[ M c^2 + \left( E_i - E_0 \right)
 * \right]^2 + E_\gamma^2}. \f]
 *
 * From this equation, the relation between the photon energy and the
 * level-energy difference is found to be:
 *
 * \f[
 *      E_\gamma = \left( E_i - E_0 \right) \left( 1 + \frac{E_i - E_0}{ 2 M
 * c^2} \right). \f]
 *
 * As can be seen, the energy of the photon needs to be higher than the
 * level-energy difference, to account for the recoil of the nucleus.
 *
 * ## Emission
 *
 * The emission process adds another recoil, which is not in the same direction
 * as the first one, in general. In addition, the emission of the gamma ray
 * during the deexcitation will take place while the nucleus is in motion, so
 * its energy will be Doppler shifted in the lab frame. Both processes will be
 * treated separately in the following.
 *
 * ### Recoil
 *
 * First, the recoil correction for the emission of a gamma ray by a nucleus at
 * rest will be derived, then the four-vector of the emitted particle will be
 * Lorentz boosted into the reference frame of the initial recoil from the
 * absorption process.
 *
 * For the decay of the nucleus in its rest frame, the four momenta before and
 * after the decay are:
 *
 * \f[
 *      p^\prime_R = \left( M c^2 + \left( E_i - E_0 \right), 0, 0, 0 \right)
 * \f]
 * \f[
 *      p^{\prime \prime}_R c = \left( \sqrt{\left[ M c^2 + \left( E_j - E_0
 * \right) \right]^2 + |\mathbf{p}^{\prime \prime}_R|^2 c^2}, \mathbf{p}^{\prime
 * \prime}_R c \right) \f] \f[ p^{\prime \prime}_{\gamma, R} c = \left(
 * E^{\prime \prime}_{\gamma, R} , \mathbf{p}^{\prime \prime}_{\gamma, R} c
 * \right) \f]
 *
 * The equations above allow for a general direction of emission and a decay to
 * any state \f$j\f$, which may be the ground state \f$j=0\f$ or any other
 * low-lying excited state
 * (\f$E_0 < E_j < E_i\f$).
 * Conservation of momentum and energy lead to require that
 *
 * \f[
 *      \mathbf{p}^{\prime \prime}_R = -\mathbf{p}^{\prime \prime}_{\gamma, R}
 * \f]
 *
 * and
 *
 * \f[
 *      M c^2 + \left( E_i - E_0 \right) = \sqrt{\left[ M c^2 + \left( E_j - E_0
 * \right) \right]^2 + E^{\prime \prime 2}_{\gamma, R}} + E^{\prime
 * \prime}_{\gamma, R} \f].
 *
 * By solving this equation for \f$E^{\prime \prime}_{\gamma, R}\f$, the
 * following expression is obtained:
 *
 * \f[
 *      E^{\prime \prime}_{\gamma, R} = \left( E_i - E_j \right) \frac{1 +
 * \frac{\left(E_i - E_0\right)^2 - \left(E_j - E_0\right)^2}{2 M c^2 \left( E_i
 * - E_j \right)}}{1 + \frac{E_i - E_0}{M c^2}}. \f]
 *
 * In this equation, it is instructive to compare the magnitude of the two
 * fractions in the numerator and denominator of the large fraction. For the
 * emission process, the gamma-ray energy is expected to be lower than the
 * level-energy difference \f$E_i - E_j\f$, i.e. the inequality
 *
 * \f[
 *      \frac{\left( E_i - E_0 \right)^2 - \left( E_j - E_0 \right)^2}{2 \left(
 * E_i - E_j \right)} < E_i - E_0 \f]
 *
 * should hold, considering that
 *
 * \f[
 *      E_0 < E_j < E_i.
 * \f]
 *
 * For the author of this text, this was not immediately obvious, but it can be
 * seen after a few manipulations. First, multiply by the denominator, and
 * expand the expressions:
 *
 * \f[
 *      \left( E_i - E_0 \right)^2 - \left( E_j - E_0 \right)^2 < 2 \left( E_i -
 * E_0 \right) \left( E_i - E_j \right) \f] \f[ E_i^2 - 2 E_i E_0 - E_j^2 + 2
 * E_j E_0 < 2 E_i^2 - 2 E_i E_0 - 2 E_j E_i + 2 E_j E_0. \f]
 *
 * The terms \f$2 E_i E_0\f$ and \f$2 E_j E_0 \f$ appear on both sides of the
 * inequality, so they can be dropped:
 *
 * \f[
 *      E_i^2 - E_j^2 < 2 E_i^2 - 2 E_j E_i.
 * \f]
 *
 * From this, the inequality
 *
 * \f[
 *      E_i^2 - 2 E_j E_i + E_j^2 = (E_i - E_j)^2 > 0
 * \f]
 *
 * can be obtained by arranging all terms to the same side.
 * Considering the assumed order of excited states, the last inequality is
 * always true. This means that the energy of the emitted gamma ray is always
 * smaller than the level energy, if the mass of the resting emitting particle
 * is finite.
 *
 * In the case of an "elastic" NRF process (\f$E_j = E_0\f$), the expression for
 * the energy of the emitted gamma ray simplifies to:
 *
 * \f[
 *      E^{\prime \prime}_{\gamma, R} = \left( E_i - E_0 \right) \frac{1 +
 * \frac{E_i - E_0}{2 M c^2}}{1 + \frac{E_i - E_0}{M c^2}}. \f]
 *
 * Although the capture of a photon by a nucleus at rest and the emission of a
 * photon by a nucleus at rest are equivalent in classical mechanics, the recoil
 * correction is different in relativistic kinematics. Using the series
 * expansion
 *
 * \f[
 *      \frac{1 + x/2}{1 + x} = \sum_{i=0}^\infty \frac{\left( -x \right)^i}{2},
 * \f]
 *
 * an expression that is more similar to the capture process can be obtained.
 * If the excitation energies are much smaller than the energy equivalent of the
 * nuclear mass
 *
 * \f[
 *      E_i - E_0 \ll M c^2,
 * \f]
 *
 * only the linear term
 *
 * \f[
 *      E^{\prime \prime}_{\gamma, R} = \left( E_i - E_0 \right) \left( 1 -
 * \frac{E_i - E_0}{2 M c^2} \right) \f]
 *
 * contributes, which is the same as the expression for \f$E_\gamma\f$. except
 * for a negative sign.
 *
 * ### Doppler Shift
 *
 * For any cascade step except the first one (regardless of whether it is a
 * capture process or spontaneous decay at rest), the emission of the gamma ray
 * occurs when the nucleus is in motion. Without loss of generality, the first
 * direction of motion was assumed to be the positive \f$z\f$ axis. A general
 * direction of motion can be obtained by an Euler-angle rotation.
 *
 * In order to get the laboratory-frame energy \f$E^{\prime \prime}_\gamma\f$
 * from the rest-frame energy \f$E^{\prime \prime}_{\gamma, R}\f$, the gamma
 * ray's four momentum needs to be Lorentz boosted along the \f$z\f$ axis. The
 * boost will modify the energy of the gamma ray in the laboratory frame, but
 * also the emission angle. The change in energy is given by the relativistic
 * Doppler effect (See, e.g., the textbook by Jackson \cite Jackson1998. In the
 * cited edition of the book, the relativistic Doppler shift is given for the
 * transformation into a moving frame. Here, the equation has been solved for
 * the energy in the laboratory frame. The author also found the dissertation of
 * C. Stahl \cite Stahl2015 and an associated publication \cite Stahl2017
 * helpful.):
 *
 * \f[
 *      E^{\prime \prime}_{\gamma} = \frac{E^{\prime \prime}_{\gamma,
 * R}}{\gamma^\prime \left[ 1 - \beta^\prime \cos \left( \theta \right) \right]}
 * \f]
 *
 * The expression for the Doppler shift of the energy uses the relativistic
 * \f$\beta\f$ and \f$\gamma\f$ factors of the recoiling nucleus after capturing
 * a photon,
 *
 *
 * \f[
 *      \beta^\prime = \frac{E_\gamma}{\sqrt{\left[ M c^2 + \left( E_i - E_0
 * \right) \right]^2 + E_\gamma^2}} \f]
 *
 * \f[
 *      \gamma^\prime = \frac{\sqrt{\left[ Mc^2 + \left( E_i - E_0 \right)
 * \right]^2 + E_\gamma^2}}{Mc^2 + \left( E_i - E_0 \right)} \f]
 *
 * and the
 * scattering angle in the laboratory frame, i.e. the angle between the
 * direction of the incoming photon and the emitted gamma ray. The conversion
 * between the scattering angle in the rest frame and the laboratory frame is
 * given by \cite Lesser2010 (The relation can also be derived from the
 * information given in Jackson \cite Jackson1998 or found in the two
 * aforementioned publications by Stahl et al., but Lesser and Cline give it
 * explicitly. Also, the author of this text wanted to point out the
 * Lesser-Cline publication since it treats the problem of angular distributions
 * in a different way: Instead of sampling the angular correlation in the rest
 * frame and transforming the four momenta, the authors transform the angular
 * correlation function. While the approach of Lesser and Cline is probably
 * computationally more efficient, it is only an approximation for low
 * velocities.):
 *
 * \f[
 *      \cos \left( \theta_R \right) = \frac{\cos \left( \theta \right) -
 * \beta^\prime}{1-\beta^\prime \cos \left( \theta \right)}. \f]
 *
 * ## Expansion for Low Velocities
 *
 * Although the small-\f$\beta^\prime\f$ expansion provides computational
 * advantages, this section has been included as a cross check for the
 * derivations in this section of the documentation.
 *
 * From the expression for \f$\beta^\prime\f$ above and the numerical examples
 * in the introduction, it can be seen that the velocity of the recoiling
 * nucleus is typically on the order of
 *
 * \f[
 *      \beta^\prime < \approx 10^{-3}.
 * \f]
 *
 * In this case, the expression for \f$E^{\prime \prime}_\gamma\f$ is often
 * expanded at first order in \f$\beta\f$. The author knows about two examples
 * in the literature in which such expansions are given explicitly for the
 * elastic NRF process, the PhD theses of Pietralla \cite Pietralla1996 and
 * Romig \cite Romig2015 (The author's own PhD thesis repeats Pietralla's
 * expression, but has an unfortunate typo.). The former is in German and not
 * publicly available. Pietralla gives the relation between \f$\left( E_i -
 * E_0\right)\f$ and \f$E^{\prime \prime}_\gamma\f$ for an elastic cascade at
 * first order in \f$E^{\prime \prime}_\gamma / M c^2\f$ as [Eq. (4.1) in \cite
 * Pietralla1996]:
 *
 * \f[
 *      E_i - E_0 = E^{\prime \prime}_\gamma \left\{ 1 + \frac{1}{2} \left[ 1 -
 * 2 \cos \left( \theta \right) \right] \frac{E^{\prime \prime}_\gamma}{M c^2}
 * \right\} + \mathcal{O} \left[ \left( \frac{E^{\prime \prime}_\gamma}{M c^2}
 * \right)^2 \right]. \f]
 *
 * Romig gives the same equation in Ref. \cite Romig2015 [Eq. (4.6)], but
 * without indicating that it is an approximation. Romig has also solved the
 * equation for \f$E_\gamma\f$ and gives [Eq. (4.5) in \cite Romig2015]:
 *
 * \f[
 *      E^{\prime \prime}_\gamma = \frac{Mc^2}{1-2 \cos \left( \theta \right)}
 * \left\{ \sqrt{1 + \frac{2 \left( E_i - E_0 \right)}{M c^2} \left[ 1 - 2 \cos
 * \left( \theta \right) \right]} - 1 \right\} \f]
 *
 * This equation has the serious flaw that it appears to be a higher-order
 * expression in \f$\left( E_i - E_0\right) / M c^2\f$ although it was derived
 * from a first-order approximation. This inconsistency causes the equation to
 * have an unphysical pole at  * \f$\theta = 60^\circ\f$. The approximative
 * character is not indicated in Ref. \cite Romig2015.
 *
 * The small-\f$\beta\f$ approximation for elastic cascades can be obtained from
 * the present expression for \f$E^{\prime \prime}_\gamma\f$. Using the relation
 * between \f$\gamma^\prime\f$ and \f$\beta^\prime\f$ yields:
 *
 * \f[
 *      E^{\prime \prime}_\gamma = E^{\prime \prime}_{\gamma, R}
 * \frac{\sqrt{1-\beta^{\prime 2}}}{1-\beta^\prime \cos \left( \theta \right)}.
 * \f]
 *
 * At first order in \f$\beta^\prime\f$, this is approximately equal to:
 *
 * \f[
 *      E^{\prime \prime}_\gamma = E^{\prime \prime}_{\gamma, R} \left[ 1 +
 * \beta^\prime \cos \left( \theta \right) \right] + \mathcal{O} \left(
 * \beta^{\prime 2} \right). \f]
 *
 * Using the relation between the gamma-ray energy in the rest frame
 * \f$E^{\prime \prime}_{\gamma,R}\f$ and the level-energy difference
 * \f$\left( E_i - E_0 \right)\f$ plus the fact that
 * \f$\beta^\prime \approx \left( E_i - E_0 \right) / M c^2\f$, and expanding
 * the equation in \f$\beta^\prime\f$ again yields:
 *
 * \f[
 *      E^{\prime \prime}_\gamma = \left( E_i - E_0 \right) \left\{ 1 -
 * \frac{1}{2} \left[ 1 - 2 \cos \left( \theta \right) \right] \frac{\left( E_i
 * - E_0 \right)}{M c^2} \right\} + \mathcal{O} \left[ \left( \frac{E_i - E_0
 * }{M c^2} \right)^2 \right]. \f]
 *
 * This expression is the same as Pietralla's with the symbols \f$\left( E_i -
 * E_0 \right)\f$ and \f$E^{\prime \prime}_\gamma\f$ swapped and an additional
 * negative sign. This makes sense, because by expanding the expressions at low
 * velocities, the classical-mechanics solution was obtained in which the
 * absorption and emission are symmetric.
 */
class FourMomentumSampler {};

} // namespace alpaca
