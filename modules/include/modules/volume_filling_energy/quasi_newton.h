#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;

namespace modules {
namespace VolumeFillingEnergy {
struct Options {
  double alpha = 1.0;
  double beta = 1.0;
  double lambda = 1.0;
  double rMax = 1.0;
  int nDirections = 32;
  bool varyingAlpha = false;
  std::string varyingAlphaType = "y"; // "y" | "circle" | "gaussian"
};
}

std::tuple<std::vector<Vector3>, std::vector<Vector3>, double>
volume_filling_energy_quasi_newton(const std::vector<Vector3> &nodes,
                      const std::vector<std::array<int, 2>> &segments,
                      const std::vector<bool> &isFixedNode,
                      const std::vector<Vector3> &meshNodes,
                      const VolumeFillingEnergy::Options &options,
                      const bool energyOnly=false);
}
