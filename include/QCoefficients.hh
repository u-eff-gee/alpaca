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
 * Q coefficients are the results of integrals involving Legendre polynomials and 
 * detector-response function that can increase the efficiency of realistic angular-correlation 
 * calculations.
 * This class description starts by introducing the problem of treating real-world detection 
 * systems.
 * After that, common methods to solve the complicated integrals that arise are introduced, and it
 * is shown how the introduction of (generalized) Q coefficients can lead to a large boost in efficiency.
 * 
 * A realistic detection system can be modeled by a response function 
 * \f$\epsilon \left( \Omega, ... \right)\f$ that depends on the solid angle 
 * \f$\Omega = \left( \theta, \varphi\right)\f$
 * (there are usually also other dependencies like the photon energy, but they will not be 
 * shown explicitly here without loss of generality because the angular correlations only depend on 
 * the solid angle.).
 * Please note that the notion of a "detection system" is very general.
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
 * Since the result of the integrations has no dependence on the solid angle, the arguments 
 * "$\left( \theta, \varphi \right)$" of the effective angular correlation should be seen as 
 * symbols indicating the position of a detector rather than variables.
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
 * respective response functions \f$\epsilon_i\f$, \f$\epsilon_j\f$.
 * In this case, the integrals over the solid angles will be restricted, and it is important to 
 * indicate whether the two photons can be distinguished.
 * For example, if the first photon is detected by detector \f$i\f$, and the second photon by 
 * detector \f$j\f$, then
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
 * If the two photons cannot be distinguished, for example because the data acquisition system 
 * cannot resolve the energy {this example is given for Eq. (4) in \cite CampvanLehn1969 , which 
 * is a practical application of Ref. \cite Rose1953. It is already formulated in terms of the Q 
 * coefficients.}, all possible contributions to an event of interest must be added up:
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
 * Depending on whether one or both detectors are polarization sensitive, even more terms may have 
 * to be summed in the equation above.
 * 
 * Although analytical and experimental methods for the determination of the effective angular 
 * correlation exist \cite Rose1953 \cite CampvanLehn1969, the complicated integral is usually
 * solved by Monte Carlo particle simulations i.e. Monte Carlo integration.
 * In the most simple case that uses isotropically sampled emission angles 
 * \f$\theta_\mathrm{rand}\f$ and \f$\varphi_\mathrm{rand}\f$, the integral in the numerator can be
 * approximated by:
 * 
 * \f[
 *       \int \int W \left( \theta_{1,2}, \varphi_{1,2} \right) \epsilon \left( \Omega_1 \right) \epsilon \left( \Omega_2 \right) \mathrm{d} \Omega_1 \mathrm{d} \Omega_2 \approx \frac{1}{N} \sum_{i = 0}^{N-1} W \left( \theta_{\mathrm{rand}, i,j}, \varphi_{\mathrm{rand}, i,j} \right) \epsilon \left( \Omega_{\mathrm{rand}, 1} \right) \epsilon \left( \Omega_{\mathrm{rand}, 2} \right).
 * \f]
 * 
 * where \f$N\f$ is the number of randomly sampled particles.
 * This is a computationally expensive method compared to the methods described below, since it 
 * requires resampling for every new angular correlation of interest.
 * 
 * Post-processing methods like the one proposed by Longland \cite Longland2010 (Sec. 4.2.2) record
 * the emission angles of the detected events from a Monte Carlo simulation with isotropic 
 * emission, and then re-weight the detected events with the angular correlation of interest.
 * Usually, the re-weighting procedure is faster than the Monte Carlo simulation, which has to be
 * performed only once in this approach.
 * 
 * Instead of performing the post-processing procedure for every angular correlation of interest,
 * it is instructive to realize that every two-step dir-dir or pol-dir angular-correlation function 
 * can be described by a finite series with a few terms that contain the angular dependence as 
 * (associated) Legendre polynomials (see, e.g., \cite Iliadis2021):
 * 
 * \f[
 *      W \left( \theta, \varphi \right) = a_0 P_0 \left[ \cos \left( \theta \right) \right] + \sum_{n = 2, 4, ...} a_n P_n \left[ \cos \left( \theta \right) \right] + b_n P_n^{|2|} \left[ \cos \left( \theta \right) \right] \cos \left( 2 \varphi \right).
 * \f]
 * 
 * The expansion coefficients \f$a_n\f$ and \f$b_n\f$ are factors that depend on the spins, 
 * parities, and mixing ratios of the angular correlation, and the time required to compute them
 * is usually negligible compared to geometrical calculations.
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
 * [using a single variable, i.e. \f$Q_n \left( \theta \right)\f$] and the pol-dir-related terms 
 * [using two variables, i.e. \f$Q_n \left( \theta, \varphi \right)\f$].
 * The dir-dir-correlation related terms can be simplified further to the product of a 
 * Legendre polynomial and the "real" angle-independent Q coefficients from the literature 
 * {\cite Rose1953 [Eq. (5), the unnumbered equation above that one, and Eq. (6)] 
 * \cite CampvanLehn1969 \cite Longland2010 }:
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
 * spin sequence and any placement of the detectors.
 * The general Q coefficients introduced here have lost the latter property because the rotational 
 * symmetry is broken by the parity measurement.
 * Nevertheless, any Q-coefficient based calculation can be expected to be a lot more efficient 
 * in cases where many different effective angular-correlation function need to be evaluated, for
 * example in a minimization procedure to determine spins, parities, and mixing ratios.
 * 
 * Due to the way unobserved intermediate transitions are treated in the formalism 
 * {see, e.g., Eqs. (22) and (23) in Ref. \cite Iliadis2021 }, i.e. by the simple insertion of 
 * another expansion coefficient, the application of Q coefficients to this case is straightforward. 
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