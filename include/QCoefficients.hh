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

    Copyright (C) 2021, 2022 Udo Friman-Gayer
*/

#pragma once

/**
 * @brief Class for generalized Q coefficients to incorporate realistic detector-response effects.
 * 
 * This class description starts by introducing the problem of treating realistic detection 
 * systems.
 * After that, common methods to solve the complicated integrals that arise are introduced, and it
 * is shown how the introduction of (generalized) Q coefficients can lead to a large boost in efficiency.
 * 
 * A realistic detection system can be modeled by a response function 
 * \f$\epsilon \left( \Omega, ... \right)\f$ that depends on the solid angle 
 * \f$\Omega = \left( \theta, \varphi\right)\f$
 * (there are usually also other dependencies like the photon energy, but they will not be 
 * shown explicitly here without loss of generality.).
 * The notion of a "detection system" is very general.
 * For example, it may also include anything that attenuates photons on their way to the actual detector.
 * 
 * Due to the existence of a response function, an effective angular correlation 
 * \f$\tilde{W} \left( \theta, \varphi \right)\f$ is observed in an experiment, which can be 
 * interpreted as a response-weighted average of the "raw" angular correlation 
 * \f$W \left( \theta, \varphi \right)\f$:
 * 
 * \f[
 *  \tilde{W} \left( \theta, \varphi \right) = 
 *      \frac{
 *          \int  \int W \left( \theta_{1,2}, \varphi_{1,2} \right) \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }{
 *          \int  \int \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }
 * \f]
 * 
 * The equation above is a generalized version of Eq. (2) from an article by M. E. Rose \cite Rose1953 
 * (not to be confused with H. J. Rose from \cite BiedenharnRose1953 and \cite RoseBrink1967).
 * The angles \f$\theta_{1,2}\f$ and \f$\varphi_{1,2}\f$ are defined by
 * 
 * \f[
 *      \begin{pmatrix}
 *          \sin \left( \theta_{1,2} \right) \cos \left( \varphi_{1,2} \right) \\
 *          \sin \left( \theta_{1,2} \right) \sin \left( \varphi_{1,2} \right) \\
 *          \cos \left( \theta_{1,2} \right)
 *      \end{pmatrix} = 
 *      \begin{pmatrix}
 *          \sin \left( \theta_{2} \right) \cos \left( \varphi_{2} \right) \\
 *          \sin \left( \theta_{2} \right) \sin \left( \varphi_{2} \right) \\
 *          \cos \left( \theta_{2} \right)
 *      \end{pmatrix} -
 *      \begin{pmatrix}
 *          \sin \left( \theta_{1} \right) \cos \left( \varphi_{1} \right) \\
 *          \sin \left( \theta_{1} \right) \sin \left( \varphi_{1} \right) \\
 *          \cos \left( \theta_{1} \right)
 *      \end{pmatrix},
 * \f]
 * 
 * where \f$\left(\theta_1, \varphi_1\right)\f$ and \f$\left(\theta_2, \varphi_2\right)\f$ are 
 * the directions of emission of the first and second photon, respectively.
 * With this choice of notation, the \f$z\f$ axis and the polarization axis are defined by the 
 * first photon, without loss of generality.
 * 
 * Usually, the detection system consists of a certain amount of detectors \f$i, j\f$ with their 
 * respective response function \f$\epsilon_{i, j}\f$.
 * In this case, the integrals over the solid angles will be restricted, and it is important to indicate whether the two photons can be distinguished, i.e., for example,
 * 
 * \f[
 *  \tilde{W} \left( \theta, \varphi \right) = 
 *      \frac{
 *          \int  \int W \left( \theta_{1,2}, \varphi_{1,2} \right) \epsilon_i \left( \Omega_1 \right) \epsilon_j \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }{
 *          \int  \int \epsilon_j \left( \Omega_1 \right) \epsilon_j \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }
 * \f]
 * 
 * or not {See, for example, Eq. (4) in \cite CampvanLehn1969 , which is a practical application 
 * of Ref. \cite Rose1953. It is already formulated in terms of the Q coefficients to be 
 * introduced here.}
 * 
 * \f[
 *  \tilde{W} \left( \theta, \varphi \right) = 
 *      \frac{
 *          \int  \int \left[ W \left( \theta_{1,2}, \varphi_{1,2} \right) \epsilon_i \left( \Omega_1 \right) \epsilon_j \left( \Omega_2 \right) + W \left( \theta_{2,1}, \varphi_{2,1} \right) \epsilon_i \left( \Omega_2 \right) \epsilon_j \left( \Omega_1 \right) \right] \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }{
 *          \int  \int \left[ \epsilon_i \left( \Omega_1 \right) \epsilon_j \left( \Omega_2 \right) + \epsilon_i \left( \Omega_2 \right) \epsilon_j \left( \Omega_1 \right) \right] \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }.
 * \f]
 * 
 * Distinguishability is even more important when polarization measurement is involved, since the
 * formalism implemented in alpaca is only capable of a single polarization measurement per
 * cascade.
 * Depending on whether one or both detectors are polarization sensitive, even more terms may have to be summed in the equation above.
 * 
 * Although analytical \cite Rose1953 and experimental \cite CampvanLehn1969 methods for the 
 * determination of the effective angular correlation exist, the complicated integral is usually
 * solved by Monte Carlo particle simulations or Monte Carlo integration.
 * In the most simple case that uses isotropically sampled emission angles 
 * \f$\theta_\mathrm{rand}\f$ and \f$\varphi_\mathrm{rand}\f$, the integral in the numerator can be
 * approximated by:
 * 
 * \f[
 *       \int \int W \left( \theta_{1,2}, \varphi_{1,2} \right) \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2 \approx \frac{1}{N} \sum_{i = 0}^{N-1} W \left( \theta_{\mathrm{rand}, i,j}, \varphi_{\mathrm{rand}, i,j} \right) \epsilon \left( \Omega_{\mathrm{rand}, 1} \right) \epsilon \left( \Omega_{\mathrm{rand}, 2} \right).
 * \f]
 * 
 * where \f$N\f$ is the number of randomly sampled particles.
 * This is probably the computationally most expensive method, since it requires resampling for 
 * every new angular correlation of interest.
 * 
 * Post-processing methods like the one proposed by Longland \cite Longland2010 (Sec. 4.2.2) record
 * the emission angles of the detected events from a Monte Carlo simulation with isotropic 
 * emission, and then re-weight the detected events with the angular correlation of interest.
 * Usually, the re-weighting procedure is faster than the Monte Carlo simulation, which has to be
 * performed only a single time in this approach.
 * 
 * Instead of performing the post-processing procedure for every angular correlation of interest,
 * it is instructive to realize that any two-step dir-dir or pol-dir angular-correlation function 
 * can be described by a finite series of terms that contain the angular dependence as (associated)
 * Legendre polynomials (see, e.g., \cite Iliadis2021):
 * 
 * \f[
 *      W \left( \theta, \varphi \right) = a_0 P_0 \left[ \cos \left( \theta \right) \right] + \sum_{n = 2, 4, ...} a_n P_n \left[ \cos \left( \theta \right) \right] + b_n P_n^{|2|} \left[ \cos \left( \theta \right) \right] \cos \left( 2 \varphi \right).
 * \f]
 * 
 * The expansion coefficients \f$a_n\f$ and \f$b_n\f$ are factors that depend on the spins and 
 * parities of the angular correlation, and they can usually be computed directly and efficiently.
 * By inserting the expression above in the general equation for the effective angular 
 * correlation, it can be seen that Monte Carlo integrations over (associated) Legendre polynomials
 * up to the required order (note that the order is limited by the involved spin quantum numbers
 * and multipolarities) can be used to represent any effective angular correlation in the given
 * detection system:
 * 
 * \f[
 *  \tilde{W} \left( \theta, \varphi \right) = 
 *      a_0 \underbrace{\frac{
 *          \int  \int P_0 \left[ \cos \left( \theta_{1,2} \right) \right] \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }{
 *          \int  \int \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }}_{Q_0 \left( \theta \right)}
 * \f]
 * \f[
 *      + \sum_{n = 2, 4, ...} a_n \underbrace{\frac{
 *          \int  \int P_n \left[ \cos \left( \theta_{1,2} \right) \right] \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }{
 *          \int  \int \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }}_{Q_n \left( \theta \right)}
 *      +b_n \underbrace{\frac{
 *          \int  \int P_n^{|2|} \left[ \cos \left( \theta_{1,2} \right) \right] \cos \left( 2 \varphi \right) \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }{
 *          \int  \int \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }}_{Q_n \left( \theta, \varphi \right)}.
 * \f]
 * 
 * In this equation, the weighted averages of the (associated) Legendre Polynomials have been 
 * denoted as generalized Q coefficients for the dir-dir-correlation-related terms 
 * [\f$Q_n \left( \theta \right)\f$] and the pol-dir-related terms 
 * [\f$Q_n \left( \theta, \varphi \right)\f$].
 * The dir-dir-correlation related terms can be simplified further to the product of a 
 * Legendre polynomial and the "real" angle-independent Q coefficients from the literature 
 * {\cite Rose1953 [Eq. (5), the unnumbered equation above that one, and Eq. (6)] \cite CampvanLehn1969 \cite Longland2010 }:
 * 
 * \f[
 *      Q_n \left( \theta \right) = \frac{
 *          \int  \int P_n \left[ \cos \left( \theta_{1,2} \right) \right] \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      }{
 *          \int  \int \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2
 *      } = Q_n P_n \left[ \cos \left( \theta \right) \right]
 * \f]
 * 
 * The real Q coefficients allow one to calculate the effective angular correlation for any
 * spin sequence and any placement of the detectors, while the general Q coefficients 
 * introduced here have lost the latter property because the rotational symmetry is broken by
 * the parity measurement.
 * Nevertheless, any Q-coefficient based calculation can be expected to be a lot more efficient 
 * in cases where many different effective angular-correlation function need to be evaluated, for
 * example in a minimization procedure to determine spins, parities, and mixing ratios.
 * 
 * Please note that the Q-coefficient formalism is possible due to the knowledge of the analytical
 * form of the dir-dir and pol-dir angular correlations.
 * Although alpaca allows to sample cascades of arbitrary length because they can be represented by
 * at most one pol-dir correlation and an arbitrary number of dir-dir correlations, this is not
 * equivalent to the knowledge of the analytical higher-order distribution function
 * \f$W \left( \theta_1, \varphi_1; \theta_2, \varphi_2, ... \right)\f$.
 */
struct QCoefficients{

};