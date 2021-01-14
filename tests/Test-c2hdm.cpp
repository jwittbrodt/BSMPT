#include "catch.hpp"

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>

const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                              /* lambda_2 = */ 0.274365,
                                              /* lambda_3 = */ 4.71019,
                                              /* lambda_4 = */ -2.23056,
                                              /* Re(lambda_5) = */ -2.43487,
                                              /* Im(lambda_5) = */ 0.124948,
                                              /* Re(m_{12}^2) = */ 2706.86,
                                              /* tan(beta) = */ 4.64487,
                                              /* Yukawa Type = */ 1};

TEST_CASE("Checking NLOVEV for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  std::vector<double> Check;
  auto sol = Minimizer::Minimize_gen_all(modelPointer,
                                         0,
                                         Check,
                                         modelPointer->get_vevTreeMin(),
                                         Minimizer::WhichMinimizerDefault);
  for (std::size_t i{0}; i < sol.size(); ++i)
  {
    auto expected = std::abs(modelPointer->get_vevTreeMin(i));
    auto res      = std::abs(sol.at(i));
    REQUIRE(std::abs(res - expected) <= 1e-4);
  }
}

TEST_CASE("Checking EWPT for C2HDM", "[c2hdm]")
{
  using namespace BSMPT;
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);
  std::vector<double> Check;
  auto EWPT = Minimizer::PTFinder_gen_all(
      modelPointer, 0, 300, Minimizer::WhichMinimizerDefault);
  const double omega_c_expected = 200.79640966130026;
  const double Tc_expected      = 145.56884765625;
  const std::vector<double> min_expected{
      0, -49.929666284908336, -194.48507400070247, 1.3351211509361505};
  REQUIRE(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS);
  REQUIRE(std::abs(omega_c_expected - EWPT.vc) / omega_c_expected <= 1e-4);
  REQUIRE(std::abs(Tc_expected - EWPT.Tc) / Tc_expected <= 1e-4);
  for (std::size_t i{0}; i < EWPT.EWMinimum.size(); ++i)
  {
    auto res      = std::abs(EWPT.EWMinimum.at(i));
    auto expected = std::abs(min_expected.at(i));
    if (expected != 0)
    {
      REQUIRE(std::abs(res - expected) / expected <= 1e-4);
    }
    else
    {
      REQUIRE(res <= 1e-4);
    }
  }
}

TEST_CASE("potential evaluation benchmark", "[c2hdm][benchmark]")
{
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      FChoose(BSMPT::ModelID::ModelIDs::C2HDM);
  modelPointer->initModel(example_point_C2HDM);

  std::vector<double> testpoint(8, 100.);

  CHECK(modelPointer->VTree(testpoint) == Approx(211159851.1352651715));
  BENCHMARK("C2HDM potential tree-level")
  {
    return modelPointer->VTree(testpoint);
  };

  CHECK(modelPointer->VEff(testpoint) == Approx(255888063.1856197417));
  BENCHMARK("C2HDM potential T=0")
  {
    return modelPointer->VEff(testpoint, 0.);
  };

  CHECK(modelPointer->VEff(testpoint, 100.) == Approx(-870825001.6903254986));
  BENCHMARK("C2HDM potential T=100.")
  {
    return modelPointer->VEff(testpoint, 100.);
  };

  CHECK(modelPointer->VTree(testpoint, 4) == Approx(-546954.4274349774));
  CHECK(modelPointer->VTree(testpoint, 6) == Approx(7271060.5550102228));
  BENCHMARK("C2HDM dV tree-level")
  {
    std::vector<double> derivative(8);
    for (std::size_t i = 0; i != 8; ++i)
    {
      derivative[i] = modelPointer->VTree(testpoint, i);
    }
    return derivative;
  };

  CHECK(modelPointer->VEff(testpoint, 0., 4) == Approx(-483348.4516666534));
  CHECK(modelPointer->VEff(testpoint, 0., 6) == Approx(6872508.8383762874));
  BENCHMARK("C2HDM dV T=0.")
  {
    std::vector<double> derivative(8);
    for (std::size_t i = 0; i != 8; ++i)
    {
      derivative[i] = modelPointer->VEff(testpoint, 0., i);
    }
    return derivative;
  };

  CHECK(modelPointer->VEff(testpoint, 100., 4) == Approx(-589695.8167427159));
  CHECK(modelPointer->VEff(testpoint, 100., 6) == Approx(6631485.4022827037));
  std::vector<double> derivative(8);
  BENCHMARK("C2HDM dV T=100.")
  {
    for (std::size_t i = 0; i != 8; ++i)
    {
      derivative[i] = modelPointer->VEff(testpoint, 100., i);
    }
  };
}
