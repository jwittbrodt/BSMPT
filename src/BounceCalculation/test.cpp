#include "VevaciousPlusPlus.hpp"

int main()
{
  std::string tunnelingStrategy       = "QuantumThenThermal";
  double survivalProbabilityThreshold = 0.01;
  int thermalIntegrationResolution    = 5;
  int temperatureAccuracy             = 7;
  int resolutionOfPathPotential       = 100;
  int pathFindingTimeout              = 10000000;
  double vacuumSeparationFraction     = 0.2;

  std::string bouncePotentialFitClass = "BubbleShootingOnPathInFieldSpace";
  std::string bouncePotentialFitArguments("<NumberShootAttemptsAllowed>"
                                          "32"
                                          "</NumberShootAttemptsAllowed>"
                                          "<RadialResolution>"
                                          "  0.05"
                                          "</RadialResolution>");

  std::string tunnelPathFinders(
      "    <TunnelPathFinders>"
      "        <PathFinder>"
      "          <ClassType>"
      "            MinuitOnPotentialOnParallelPlanes"
      "          </ClassType>"
      "          <ConstructorArguments>"
      "            <NumberOfPathSegments>"
      "              50"
      "            </NumberOfPathSegments>"
      "            <MinuitStrategy>"
      "              1"
      "            </MinuitStrategy>"
      "            <MinuitTolerance>"
      "              0.5"
      "            </MinuitTolerance>"
      "          </ConstructorArguments>"
      "        </PathFinder>"
      "        <PathFinder>"
      "          <ClassType>"
      "            MinuitOnPotentialPerpendicularToPath"
      "          </ClassType>"
      "          <ConstructorArguments>"
      "            <NumberOfPathSegments>"
      "              100"
      "            </NumberOfPathSegments>"
      "            <NumberOfAllowedWorsenings>"
      "              1"
      "            </NumberOfAllowedWorsenings>"
      "            <ConvergenceThresholdFraction>"
      "              0.05"
      "            </ConvergenceThresholdFraction>"
      "            <MinuitDampingFraction>"
      "              0.75"
      "            </MinuitDampingFraction>"
      "            <NeighborDisplacementWeights>"
      "              0.5"
      "              0.25"
      "            </NeighborDisplacementWeights>"
      "            <MinuitStrategy>"
      "              1"
      "            </MinuitStrategy>"
      "            <MinuitTolerance>"
      "              0.5"
      "            </MinuitTolerance>"
      "          </ConstructorArguments>"
      "        </PathFinder>"
      "      </TunnelPathFinders>");

  auto pathFinders =
      VevaciousPlusPlus::VevaciousPlusPlus::CreateBouncePathFinders(
          tunnelPathFinders);

  auto bounceActionCalculator =
      VevaciousPlusPlus::VevaciousPlusPlus::CreateBounceActionCalculator(
          bouncePotentialFitClass, bouncePotentialFitArguments);

  std::unique_ptr<VevaciousPlusPlus::BounceAlongPathWithThreshold>
      bouncerAlongPathWithThreshold = VevaciousPlusPlus::Utils::make_unique<
          VevaciousPlusPlus::BounceAlongPathWithThreshold>(
          std::move(pathFinders),
          std::move(bounceActionCalculator),
          VevaciousPlusPlus::VevaciousPlusPlus::InterpretTunnelingStrategy(
              tunnelingStrategy),
          survivalProbabilityThreshold,
          thermalIntegrationResolution,
          temperatureAccuracy,
          resolutionOfPathPotential,
          pathFindingTimeout,
          vacuumSeparationFraction);
  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}
