/*
 * ThermalFunctions.h
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 */

#ifndef INCLUDE_BSMPT_THERMALFUNCTIONS_THERMALFUNCTIONS_H_
#define INCLUDE_BSMPT_THERMALFUNCTIONS_THERMALFUNCTIONS_H_
namespace BSMPT
{
namespace ThermalFunctions
{

/**
 * Value for m^2/T^2 at which the interpolated is changed between the two
 * polynomials, see x_+^2 in Eq. (2.42) in the manual
 */
constexpr double C_FermionTheta = 2.216079120;
/**
 * constant shift to the polynomial to make the expansions continuous, see
 * delta_+ in Eq. (2.42) in the manual
 */
constexpr double C_FermionShift = -0.01560316619;

/**
 * Value for m^2/T^2 at which the interpolated is changed between the two
 * polynomials, see x_-^2 in Eq. (2.43) in the manual
 */
constexpr double C_BosonTheta = 9.469230596;
/**
 * constant shift to the polynomial to make the expansions continuous, see
 * delta_- in Eq. (2.43) in the manual
 */
constexpr double C_BosonShift = 0.0063108787;

/**
 * @brief C_euler_gamma Euler gamma constant
 */
constexpr double C_euler_gamma = 0.5772156649015328606065;

/**
 * Integrand of the thermic integral for the bosons \f$ J_-(x) =
 * \int\limits_{0}^{\infty} \,\mathrm{d}k \, k^2 \log\left[ 1 - \exp\left(
 * -\sqrt{k^2+x} \right) \right] \f$
 * @param x The ratio m^2/T^2
 * @param k The integration variable
 * @param diff Returns the integrand of J_- for diff = 0 and for the dJ_-/dx for
 * diff = 1
 */
double JbosonIntegrand(double x, double k, int diff = 0);

/**
 * Numerical integration of the thermical integral for the bosons \f$ J_-(x) =
 * \int\limits_{0}^{\infty} \,\mathrm{d}k \, k^2 \log\left[ 1 - \exp\left(
 * -\sqrt{k^2+x} \right) \right] \f$
 * @param x The ratio m^2/T^2
 * @param diff Returns the numerical integration of J_- for diff = 0 and J_-/dx
 * for diff = 1
 */
double JbosonNumericalIntegration(double x, int diff = 0);

/**
 * Taylor expansion of J_- for small x=m^2/T^2, \f$ J_{_,s}(x,n) =
 * -\frac{\pi^4}{45} + \frac{\pi}{12} x - \frac{\pi}{6} x^{3/2} -
 * \frac{x^2}{32}(\log x - c_-) + \pi^2 x \sum\limits_{l=2}^{n} \left( -
 * \frac{x}{4\pi^2} \right)^l \frac{(2l-3)!!\zeta(2l-1)}{(2l)!! (l+1)}  \f$  ,
 * c.f. Eq. (2.38) in the manual
 * @param x The ratio m^2/T^2
 * @param n The order of the taylor expansion
 * @param diff Returns the expansion for diff = 0 and its derivative for diff =
 * 1
 */
double JbosonInterpolatedLow(double x, int n, int diff = 0);

/**
 * Using linear interpolation with data points to interpolate the thermal
 * integral for bosons for x=m^2/T^2 < 0
 * @param x The ratio m^2/T^2
 */
double JbosonInterpolatedNegative(double x);

/**
 * Puts together the separate interpolations for J_-
 * @param x The ratio m^2/T^2
 * @param diff Returns the interpolation of J_- for diff = 0 and for dJ_-/dx for
 * diff = 1
 */
double JbosonInterpolated(double x, int diff = 0);

/**
 * Integrand of the thermic integral for the fermions \f$ J_+(x) =
 * \int\limits_{0}^{\infty} \,\mathrm{d}k \, k^2 \log\left[ 1 + \exp\left(
 * -\sqrt{k^2+x} \right) \right] \f$
 * @param x The ratio m^2/T^2
 * @param k The integration variable
 * @param diff Returns the integrand of J_+ for diff = 0 and for the dJ_+/dx for
 * diff = 1
 */
double JfermionIntegrand(double x, double k, int diff = 0);

/**
 * Numerical integration of the thermical integral for the fermions \f$ J_+(x) =
 * \int\limits_{0}^{\infty} \,\mathrm{d}k \, k^2 \log\left[ 1 + \exp\left(
 * -\sqrt{k^2+x} \right) \right] \f$
 * @param x The ratio m^2/T^2
 * @param diff Returns the integrand of J_+ for diff = 0 and for the dJ_+/dx for
 * diff = 1
 */
double JfermionNumericalIntegration(double x, int diff = 0);

/**
 * Taylor expansion of J_+ for small x=m^2/T^2, \f$ J_{+,s}(x,n) =
 * -\frac{7\pi^4}{360} + \frac{\pi^2}{24} x + \frac{x^2}{32}\left( \log x - c_+
 * \right) - \pi^2 x \sum\limits_{l=2}^{n} \left( - \frac{x}{4\pi^2} \right)^l
 * \frac{(2l-3)!!\zeta(2l-1)}{(2l)!! (l+1)}  \left( 2^{2l-1} - 1\right)  \f$  ,
 * see J_{+,s} in Eq. (2.37) in the manual
 * @param x The ratio m^2/T^2
 * @param n The order of the taylor expansion
 * @param diff Returns the expansion for diff = 0 and dJ_+/dx for diff = 1
 */
double JfermionInterpolatedLow(double x, int n, int diff = 0);

/**
 * Puts together the separate interpolations for J_+, see Eq. (2.44) in the
 * manual
 * @param x The ratio m^2/T^2
 * @param diff Returns the interpolation of J_+ for diff = 0 and dJ_+/dx for
 * diff = 1
 */
double JfermionInterpolated(double x, int diff = 0);

/**
 * Expansion for large x = m^2/T^2 of the thermal integrals, \f$ J_{\pm,l}(x,n)
 * = -\exp\left(x^{1/2}\right) \left( \frac{\pi}{2} x^{3/2} \right)^{1/2}
 * \sum\limits_{l=0}^{n} \frac{1}{2^l l!} \frac{\Gamma(5/2+l)}{\Gamma(5/2-l)}
 * x^{-l/2} \f$, cf. Eq. (2.41) in the manual
 * @param x The ratio m^2/T^2
 * @param n The order of the expansion
 * @param diff Returns the expansion for diff = 0 and for its derivative for
 * diff = 1
 */
double JInterpolatedHigh(double x, int n, int diff = 0);

} // namespace ThermalFunctions
} // namespace BSMPT

#endif /* INCLUDE_BSMPT_THERMALFUNCTIONS_THERMALFUNCTIONS_H_ */
