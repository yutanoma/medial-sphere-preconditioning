#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;

namespace modules {
std::tuple<std::vector<std::vector<Vector3>>, std::vector<std::vector<int>>>
medial_axis_3d(const std::vector<Vector3> &nodes,
               const std::vector<std::array<int, 2>> &segments,
               const std::vector<Vector3> &meshNodes,
               const double maxRadius,
               const int nDirections);
}