#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;

namespace modules {
// remesh the curve such that all the edge lengths are between h and 2h
std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>, std::vector<bool>>
remesh_curve(const std::vector<Vector3> &nodes,
             const std::vector<std::array<int, 2>> &segments,
             const std::vector<bool> &isFixed,
             double h);
}
