/*
 * ThermalFunctions.cpp
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

#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/models/SMparam.h>

#include <array>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

namespace BSMPT
{
namespace ThermalFunctions
{

namespace
{
constexpr double piSquared     = M_PI * M_PI;
constexpr double piToTheFourth = M_PI * M_PI * M_PI * M_PI;
} // namespace

double JbosonIntegrand(double x, double k, int diff)
{
  auto root = std::sqrt(std::complex<double>(k * k + x));
  if (diff == 0)
  {
    return std::pow(k, 2) * std::log(1. - std::exp(-root)).real();
  }
  else if (diff == 1)
  {
    return std::pow(k, 2) *
           (std::exp(-root) / (2. * root * (1. - std::exp(-root)))).real();
  }
  return 0.;
}

double JfermionIntegrand(double x, double k, int diff)
{
  auto root = std::sqrt(std::complex<double>(k * k + x));
  if (diff == 0)
  {
    return std::pow(k, 2) * std::log(1. + std::exp(-root)).real();
  }
  else if (diff == 1)
  {
    return -std::pow(k, 2) *
           (exp(-root) / (2. * root * (1. + exp(-root)))).real();
  }
  return 0.;
}

double JbosonNumericalIntegration(double x, int diff)
{
  const std::size_t workspaceSize = 1000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(workspaceSize);
  gsl_function F;
  struct GSLTemp
  {
    double x;
    int diff;
  };
  GSLTemp dummy;
  dummy.diff = diff;
  dummy.x    = x;
  F.function = [](double k, void *p) {
    struct GSLTemp *params = static_cast<GSLTemp *>(p);
    return JbosonIntegrand(params->x, k, params->diff);
  };
  F.params = &dummy;
  double error, result;
  double upperLimit = 30;
  gsl_integration_qags(
      &F, 0, upperLimit, 0, 1e-7, workspaceSize, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double JfermionNumericalIntegration(double x, int diff)
{
  const std::size_t workspaceSize = 1000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(workspaceSize);
  gsl_function F;
  struct GSLTemp
  {
    double x;
    int diff;
  };
  GSLTemp dummy;
  dummy.diff = diff;
  dummy.x    = x;
  F.function = [](double k, void *p) {
    struct GSLTemp *params = static_cast<GSLTemp *>(p);
    return JfermionIntegrand(params->x, k, params->diff);
  };
  F.params = &dummy;
  double error, result;
  double upperLimit = 30;
  gsl_integration_qags(
      &F, 0, upperLimit, 0, 1e-7, workspaceSize, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double JfermionInterpolatedLow(double x, int n, int diff)
{
  // (2*l-3)!! * zeta(2*l-1) / ((2*l)!! * (l+1)) * (pow(2, 2*l-1) - 1)
  // first two entries are not needed and set to zero
  static constexpr std::array<double, 20> Kl{0.,
                                             0.,
                                             3.50599930088215050e-01,
                                             5.02261881397569887e-01,
                                             1.00047154865237675e+00,
                                             2.33345313876230254e+00,
                                             6.00003399630653966e+00,
                                             1.65000103629144554e+01,
                                             4.76666699902108491e+01,
                                             1.43000001107512020e+02,
                                             4.42000000380316351e+02,
                                             1.39966666680047683e+03,
                                             4.52200000004803405e+03,
                                             1.48580000000175387e+04,
                                             4.95266666666731617e+04,
                                             1.67152500000002445e+05,
                                             5.70285000000000931e+05,
                                             1.96431500000000023e+06,
                                             6.82340999999999907e+06,
                                             2.38819350000000000e+07};

  if (x == 0 and diff == 0)
  {
    return -7 * piToTheFourth / 360.0;
  }
  else if (x == 0 and diff == 1)
  {
    return piSquared / 24;
  }
  using std::log;
  using std::pow;
  double res = 0;
  static const double cf =
      1.5 + 2 * log(4 * M_PI) - 2 * C_euler_gamma - 2 * log(4);
  if (diff == 0)
  {
    res = -7 * piToTheFourth / 360.0 + piSquared / 24 * x +
          1 / 32.0 * pow(x, 2) * (log(x) - cf);
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      sum += pow(-x / (4 * piSquared), l) * Kl[l];
    }
    res += -piSquared * x * sum;
  }
  else if (diff == 1)
  {
    res        = piSquared / 24.0 + x * (-6 * cf + 3) / 96 + x * log(x) / 16;
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      sum += -Kl[l] * pow(-x / 4.0, l) * (l + 1) * pow(M_PI, 2 - 2 * l);
    }
    res += sum;
  }

  return res;
}

double JbosonInterpolatedLow(double x, int n, int diff)
{
  // (2*l-3)!! * zeta(2*l-1) / ((2*l)!! * (l+1))
  // first two entries are not needed and set to zero
  static constexpr std::array<double, 20> Kl{0.,
                                             0.,
                                             5.00857042983164358e-02,
                                             1.62019961741151561e-02,
                                             7.87772872954627286e-03,
                                             4.56644449855636483e-03,
                                             2.93113531817613092e-03,
                                             2.01440732058533233e-03,
                                             1.45471571978548087e-03,
                                             1.09101175017747658e-03,
                                             8.43049704418222007e-04,
                                             6.67413394076285765e-04,
                                             5.39064471615851638e-04,
                                             4.42802919233454986e-04,
                                             3.69002424446311487e-04,
                                             3.11345793886758830e-04,
                                             2.65559647355955364e-04,
                                             2.28676362921096933e-04,
                                             1.98587367782560958e-04,
                                             1.73763946805947944e-04};
  if (x == 0 and diff == 0)
  {
    return -piToTheFourth / 45.0;
  }
  else if (x == 0 and diff == 1)
  {
    return piSquared / 12.0;
  }
  using std::log;
  using std::pow;
  using std::sqrt;
  static const double cb = 1.5 + 2 * std::log(4 * M_PI) - 2 * C_euler_gamma;
  double res             = 0;
  if (diff == 0)
  {
    res = -piToTheFourth / 45.0;
    res += piSquared * x / 12.0;
    res += -M_PI * pow(x, 1.5) / 6;
    res += -pow(x, 2) * (log(x) - cb) / 32.0;
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      sum += pow(-x / (4 * piSquared), l) * Kl[l];
    }
    res += piSquared * x * sum;
  }
  else if (diff == 1)
  {
    res = piSquared;
    res += x * (6 * cb - 3) / 96.0;
    res += -x * log(x) / 16.0;
    res += -M_PI * sqrt(x) / 4.0;
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      sum += Kl[l] * pow(-x / 4.0, l) * (l + 1) * pow(M_PI, 2 - 2 * l);
    }
    res += sum;
  }
  return res;
}

double JInterpolatedHigh(double x, int n, int diff)
{
  // 1 / (pow(2, l) * l!) * gamma(2.5 + l) / gamma(2.5 - l)
  static constexpr std::array<double, 20> Kl{1.,
                                             1.875,
                                             8.203125e-01,
                                             -3.076171875e-01,
                                             3.17230224609375278e-01,
                                             -5.15499114990234708e-01,
                                             1.12765431404113836e+00,
                                             -3.08091267943382130e+00,
                                             1.00611054687760770e+01,
                                             -3.81483582357759801e+01,
                                             1.64514794891783822e+02,
                                             -7.94531679875092095e+02,
                                             4.24577866433252620e+03,
                                             -2.48623000632548901e+04,
                                             1.58275178081256512e+05,
                                             -1.08814184930863883e+06,
                                             8.03354724684893526e+06,
                                             -6.33823249696243107e+07,
                                             5.32147436724137306e+08,
                                             -4.73681238084051323e+09};

  using std::exp;
  using std::pow;
  using std::sqrt;

  double res = 0;
  if (diff == 0)
  {
    double sum = 0;
    for (int l = 0; l <= n; l++)
    {
      sum += Kl[l] * pow(x, -l / 2.0);
    }
    res = -exp(-sqrt(x)) * sqrt(M_PI / 2 * pow(x, 1.5)) * sum;
  }
  else if (diff == 1)
  {
    double sum = 0;
    for (int l = 0; l <= n; l++)
    {
      sum += Kl[l] * pow(x, (1 - l) / 2) * (2 * l + 2 * sqrt(x) - 3);
    }
    res = exp(-sqrt(x)) * sqrt(2 * M_PI) / (8 * pow(x, 3.0 / 4.0)) * sum;
  }
  return res;
}

double JfermionInterpolated(double x, int diff)
{
  double res = 0;
  if (x >= C_FermionTheta)
  {
    res = -JInterpolatedHigh(x, 3, diff);
  }
  else
  {
    res = -JfermionInterpolatedLow(x, 4, diff);
    if (diff == 0) res += -C_FermionShift;
  }
  return res;
}

double JbosonInterpolated(double x, int diff)
{
  double res = 0;
  if (x >= C_BosonTheta)
  {
    res = JInterpolatedHigh(x, 3, diff);
  }
  else if (x >= 0)
  {
    res = JbosonInterpolatedLow(x, 3, diff);
    if (diff == 0) res += C_BosonShift;
  }
  else if (x < 0 and diff == 0)
  {
    res = JbosonInterpolatedNegative(x);
  }
  return res;
}

double JbosonInterpolatedNegative(double x)
{
  if (x >= 0) return 0;
  double PotVal = 0;
  double xprev, fprev, xnext, fnext;
  if (-x >= C_NegLine - 2)
  {
    xprev = NegLinearInt[C_NegLine - 2][0];
    xnext = NegLinearInt[C_NegLine - 1][0];
    fprev = NegLinearInt[C_NegLine - 2][1];
    fnext = NegLinearInt[C_NegLine - 1][1];
  }
  else
  {
    std::size_t pos = (std::size_t(-x));
    xprev           = NegLinearInt[pos][0];
    fprev           = NegLinearInt[pos][1];
    xnext           = NegLinearInt[pos + 1][0];
    fnext           = NegLinearInt[pos + 1][1];
  }

  PotVal = (fnext - fprev) / (xnext - xprev) * (x - xprev) + fprev;
  return PotVal + C_BosonShift;
}

} // namespace ThermalFunctions
} // namespace BSMPT
