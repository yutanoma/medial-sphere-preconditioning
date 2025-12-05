#include "modules/line_search.h"

#include "modules/remesh_curve.h"
#include "modules/medial_axis_3d.h"
#include <igl/copyleft/cgal/is_self_intersecting.h>
#include "modules/volume_filling_energy/quasi_newton.h"
#include "polyscope/surface_mesh.h"
namespace modules {
bool ccd(const std::vector<Vector3> &nodes,
         const std::vector<Vector3> &newNodes,
        const std::vector<std::array<int, 2>> &segments,
        const std::vector<bool> &isFixedNode) {
  Eigen::MatrixXd V(nodes.size() * 2, 3);
  Eigen::MatrixXi F(segments.size() * 2, 3);

  std::vector<std::vector<int>> v2e(nodes.size());
  for (int i = 0; i < segments.size(); i++) {
    v2e[segments[i][0]].push_back(i);
    v2e[segments[i][1]].push_back(i);
  }

  std::vector<std::vector<int>> e2e(segments.size(), std::vector<int>());
  for (int i = 0; i < segments.size(); i++) {
    for (int j = 0; j < 2; j++) {
      int v = segments[i][j];
      for (int k = 0; k < v2e[v].size(); k++) {
        int e = v2e[v][k];
        if (e != i) {
          e2e[i].push_back(e);
        }
      }
    }
  }

  for (int i = 0; i < nodes.size(); i++) {
    V.row(i) << nodes[i].x, nodes[i].y, nodes[i].z;
    V.row(i + nodes.size()) << newNodes[i].x, newNodes[i].y, newNodes[i].z;
  }

  for (int i = 0; i < segments.size(); i++) {
    F.row(i) << segments[i][0], segments[i][1], segments[i][0] + nodes.size();
    F.row(i + segments.size()) << segments[i][1], segments[i][1] + nodes.size(),
        segments[i][0] + nodes.size();
  }

  auto psm = polyscope::registerSurfaceMesh("swept mesh", V, F)->setEnabled(false);

  // return false;

  // orient the faces
  // Eigen::MatrixXi FF;
  // Eigen::VectorXi C;
  // igl::bfs_orient(F, FF, C);
  // F = FF;

  std::cout << "detecting self intersections" << std::endl;

  igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
  params.detect_only = true;
  Eigen::MatrixXi IF;
  Eigen::VectorXi J,IM;
  {
    Eigen::MatrixXd tempV;
    Eigen::MatrixXi tempF;
    igl::copyleft::cgal::remesh_self_intersections(V, F, params, tempV, tempF,
                                                   IF, J, IM);
  }
  std::cout << "done self intersection detection" << std::endl;
  // std::cout << IF << std::endl;

  // Eigen::VectorXd selfintersecting = Eigen::VectorXd::Zero(F.rows());
  // for (int i = 0; i < IF.rows(); i++) {
  //   for (int j = 0; j < IF.cols(); j++) {
  //     selfintersecting(IF(i, j)) = 1;
  //   }
  // }
  // std::cout << selfintersecting << std::endl;
  // polyscope::registerSurfaceMesh("swept mesh", V, F)
  //     ->addFaceScalarQuantity("self intersecting", selfintersecting);

  // std::cout << "IF: " << IF << std::endl;

  // if IF is only happening on incident faces, it's okay
  for (int i = 0; i < IF.rows(); i++) {
    int f1 = IF(i, 0);
    int f2 = IF(i, 1);
    int e1 = f1 >= segments.size() ? f1 - segments.size() : f1;
    int e2 = f2 >= segments.size() ? f2 - segments.size() : f2;

    // std::cout << "f1: " << f1 << " f2: " << f2 << std::endl;
    // std::cout << "e1: " << e1 << " e2: " << e2 << std::endl;
    // std::cout << "e2e[e1]: " << e2e[e1].size() << std::endl;
    // std::cout << "e2e[e2]: " << e2e[e2].size() << std::endl;

    bool found = false;

    for (int j = 0; j < e2e[e1].size(); j++) {
      if (e2e[e1][j] == e2) {
        found = true;
        break;
      }
    }

    if (!found) {
      return true;
    }
  }

  std::cout << "no self intersections" << std::endl;

  return false;
}

std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>,
           std::vector<bool>>
line_search(const std::vector<Vector3> &nodes,
            const std::vector<std::array<int, 2>> &segments,
            const std::vector<bool> &isFixedNode,
            const std::vector<Vector3> &meshNodes,
            const VolumeFillingEnergy::Options &options,
            const std::vector<Vector3> &d,
            const double h,
            const double stepsize) {
  auto newNodes = nodes;
  auto newSegments = segments;

  double alpha = stepsize;

  // double totalLength = 0;
  // for (int i = 0; i < segments.size(); i++) {
  //   totalLength += (nodes[segments[i][0]] - nodes[segments[i][1]]).norm();
  // }

  // auto [_d, _g, _f] = modules::volume_filling_energy_quasi_newton(
  //     nodes, segments, isFixedNode, meshNodes, options, true);
  // double oldEnergy = _f / totalLength;

  for (int itr = 0; itr < 40; itr++) {
    auto newNodes = nodes;
    auto newSegments = segments;
    auto newIsFixed = isFixedNode;

    for (int i = 0; i < nodes.size(); i++) {
      if (isFixedNode[i]) {
        continue;
      }

      newNodes[i] = nodes[i] + alpha * d[i];
    }

    if (ccd(nodes, newNodes, segments, isFixedNode)) {
      alpha *= 0.5;
      continue;
    }

    // auto [d, g, f] = modules::volume_filling_energy_quasi_newton(
    //     newNodes, segments, isFixedNode, meshNodes, options, true);

    // double totalLength = 0;
    // for (int i = 0; i < segments.size(); i++) {
    //   totalLength += (newNodes[segments[i][0]] - newNodes[segments[i][1]]).norm();
    // }
    // double newEnergy = f / totalLength;

    // std::cout << "newEnergy: " << newEnergy << " oldEnergy: " << oldEnergy << std::endl;

    // if (newEnergy > oldEnergy) {
    //   alpha *= 0.5;
    //   continue;
    // }

    // std::tie(newNodes, newSegments, newIsFixed) =
    //     remesh_curve(newNodes, segments, isFixedNode, h);

    return {newNodes, newSegments, newIsFixed};
  }

  return {nodes, segments, isFixedNode};
}
}

