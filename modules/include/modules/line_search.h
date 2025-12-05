#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "modules/volume_filling_energy/quasi_newton.h"
using namespace geometrycentral;

namespace modules {
std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>,
           std::vector<bool>>
line_search(const std::vector<Vector3> &nodes,
            const std::vector<std::array<int, 2>> &segments,
            const std::vector<bool> &isFixedNode,
            const std::vector<Vector3> &meshNodes,
            const VolumeFillingEnergy::Options &options,
            const std::vector<Vector3> &d,
            const double h,
            const double stepsize);
}

