#include "geometrycentral/numerical/linear_algebra_types.h"
#include "modules/volume_filling_energy/newton_test.h"
#include "modules/medial_axis_3d.h"

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <cstddef>
#include <map>

#include <chrono>

#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

namespace modules {
std::tuple<std::vector<Vector3>, std::vector<Vector3>, double>
volume_filling_energy_newton_test(const std::vector<Vector3> &nodes,
                                const std::vector<std::array<int, 2>> &segments,
                                const std::vector<bool> &isFixedNode,
                                const std::vector<Vector3> &meshNodes,
                                const VolumeFillingEnergy::Options &options) {
    std::cout << "ndirections: " << options.nDirections << std::endl;
  
  auto start = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<int>> node2Segments(nodes.size(), std::vector<int>{});
  for (int i = 0; i < segments.size(); i++) {
    node2Segments[segments[i][0]].emplace_back(i);
    node2Segments[segments[i][1]].emplace_back(i);
  }

  std::vector<std::vector<int>> node2Nodes(nodes.size(), std::vector<int>{});
  for (int i = 0; i < segments.size(); i++) {
    node2Nodes[segments[i][0]].emplace_back(segments[i][1]);
    node2Nodes[segments[i][1]].emplace_back(segments[i][0]);
  }

  std::map<int, int> activeNode2NodeIdx, node2ActiveNodeIdx;
  for (int i = 0; i < nodes.size(); i++) {
    if (!isFixedNode[i]) {
      activeNode2NodeIdx[activeNode2NodeIdx.size()] = i;
      node2ActiveNodeIdx[i] = node2ActiveNodeIdx.size();
    }
  }
  int numActiveNodes = activeNode2NodeIdx.size();

  std::vector<double> segmentLengths(segments.size(), 0.0);
  double totalLength = 0.0;
  for (size_t i = 0; i < segments.size(); i++) {
    segmentLengths[i] = (nodes[segments[i][0]] - nodes[segments[i][1]]).norm();
    totalLength += segmentLengths[i];
  }

  std::vector<double> nodeWeights(nodes.size(), 0.0);
  for (int i = 0; i < nodes.size(); i++) {
    for (int j = 0; j < node2Segments[i].size(); j++) {
      int segmentId = node2Segments[i][j];
      double length = segmentLengths[segmentId];
      nodeWeights[i] += length / 2;
    }
  }

  // get the segments with 2 active nodes and 1 active node
  std::vector<std::pair<int, int>> segmentsWith2ActiveNodes = {};
  std::vector<int> activeTwoSegment2SegmentIdx = {};
  // first node is active while the second node is fixed
  std::vector<std::pair<int, int>> segmentsWith1ActiveNode = {};
  std::vector<int> activeOneSegment2SegmentIdx = {};

  for (size_t i = 0; i < segments.size(); i++) {
    int v0 = segments[i][0], v1 = segments[i][1];

    if (isFixedNode[v0] && isFixedNode[v1]) {
      continue;
    } else if (!isFixedNode[v0] && !isFixedNode[v1]) {
      int newId = segmentsWith2ActiveNodes.size();
      segmentsWith2ActiveNodes.emplace_back(std::pair<int, int> {v0, v1});
      activeTwoSegment2SegmentIdx.emplace_back(i);
    } else {
      int activeNode = isFixedNode[v0] ? v1 : v0;
      int fixedNode = isFixedNode[v0] ? v0 : v1;
      segmentsWith1ActiveNode.emplace_back(std::pair<int, int> {activeNode, fixedNode});
      activeOneSegment2SegmentIdx.emplace_back(i);
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "  preparing O(n) data: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // get the medial axis positions
  auto [medialAxisPositions, otherNodes] = modules::medial_axis_3d(
      nodes, segments, meshNodes, options.rMax, options.nDirections);

  std::vector<std::pair<int, int>> active2MedialAxisPairs;
  std::vector<Vector3> active2MedialAxisPositions;
  std::vector<std::pair<int, Vector3>> active1MedialAxes;
  for (int i = 0; i < medialAxisPositions.size(); i++) {
    for (int j = 0; j < medialAxisPositions[i].size(); j++) {
      int otherNode = otherNodes[i][j];
      if (otherNode != -1) {
        if (i != 0) {
          continue;
        }
        std::cout << "i: " << i << ", otherNode: " << otherNode << std::endl;
        active2MedialAxisPairs.emplace_back(std::pair<int, int> {i, otherNode});
        active2MedialAxisPositions.emplace_back(medialAxisPositions[i][j]);
      } else {
        // active1MedialAxes.emplace_back(std::pair<int, Vector3> {i, medialAxisPositions[i][j]});
      }
    }
  }

  std::cout << "active1MedialAxes.size(): " << active1MedialAxes.size() << std::endl;

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  sampling medial axis on S^n: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  start = std::chrono::high_resolution_clock::now();

  {
    std::vector<Vector3> medialAxisPositionsFlat;
    std::vector<Vector3> rayNodes;
    std::vector<std::vector<int>> raySegments;
    for (int i = 0; i < medialAxisPositions.size(); i++) {
      medialAxisPositionsFlat.insert(medialAxisPositionsFlat.end(),
                                     medialAxisPositions[i].begin(),
                                     medialAxisPositions[i].end());
      for (int j = 0; j < medialAxisPositions[i].size(); j++) {
        rayNodes.emplace_back(medialAxisPositions[i][j]);
        rayNodes.emplace_back(nodes[i]);
        raySegments.emplace_back(std::vector<int>{int(rayNodes.size() - 2),
                                                  int(rayNodes.size() - 1)});
        // std::cout << medialAxisPositions[i][j] << " ";
      }
      // std::cout << std::endl;
    }

    polyscope::registerPointCloud("medialAxisPositions",
                                  medialAxisPositionsFlat)->setEnabled(false);
    polyscope::registerCurveNetwork("raySegments", rayNodes, raySegments)
        ->setEnabled(false);

    polyscope::registerPointCloud("active2MedialAxes",
                                  active2MedialAxisPositions)
        ->setEnabled(false);

    std::vector<std::vector<int>> maPairs(active2MedialAxisPairs.size());
    for (int i = 0; i < active2MedialAxisPairs.size(); i++) {
      maPairs[i] = {active2MedialAxisPairs[i].first,
                    active2MedialAxisPairs[i].second};
    }

    polyscope::registerCurveNetwork("maPairs", nodes, maPairs)
        ->setEnabled(false);

    std::vector<std::vector<int>> ma1Pairs(active1MedialAxes.size());
    std::vector<Vector3> ma1Nodes = {};
    std::vector<Vector3> ma1NodesFlat = {};
    for (int i = 0; i < active1MedialAxes.size(); i++) {
      auto node = active1MedialAxes[i].second;
      auto otherNode = nodes[active1MedialAxes[i].first];
      ma1Nodes.emplace_back(node);
      ma1Nodes.emplace_back(otherNode);
      ma1NodesFlat.emplace_back(node);
      ma1Pairs[i] = {int(ma1Nodes.size() - 2), int(ma1Nodes.size() - 1)};
    }

    polyscope::registerCurveNetwork("ma1Pairs", ma1Nodes, ma1Pairs)
        ->setEnabled(false);

    polyscope::registerPointCloud("ma1NodesFlat", ma1NodesFlat)->setEnabled(false);
  }

  // Eigen::SparseMatrix<double> bilaplacian = L * L;
  
  auto func = TinyAD::scalar_function<3>(TinyAD::range(numActiveNodes));

  // step 1: add the curve length energy
  // 2 active nodes
  func.add_elements<2>(TinyAD::range(segmentsWith2ActiveNodes.size()),
                       [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;
    int segmentId = activeTwoSegment2SegmentIdx[e_id];

    int _v0 = segmentsWith2ActiveNodes[e_id].first;
    int _v1 = segmentsWith2ActiveNodes[e_id].second;

    Eigen::Vector3<T> v0 = element.variables(node2ActiveNodeIdx[_v0]);
    Eigen::Vector3<T> v1 = element.variables(node2ActiveNodeIdx[_v1]);
  
    T l = (v0 - v1).norm();
    double l_0 = segmentLengths[segmentId];

    // std::cout << "l: " << l << ", l_0: " << l_0 << std::endl;

    return options.beta * l * l / l_0;
  });

  func.add_elements<1>(TinyAD::range(segmentsWith1ActiveNode.size()),
                       [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;
    int segmentId = activeOneSegment2SegmentIdx[e_id];

    int _v0 = segmentsWith1ActiveNode[e_id].first;
    int _v1 = segmentsWith1ActiveNode[e_id].second;

    Eigen::Vector3<T> v0 = element.variables(node2ActiveNodeIdx[_v0]);
    Eigen::Vector3d v1(
      nodes[_v1].x,
      nodes[_v1].y,
      nodes[_v1].z
    );

    T l = (v0 - v1).norm();
    double l_0 = segmentLengths[segmentId];

    // std::cout << "l: " << l << ", l_0: " << l_0 << std::endl;

    return options.beta * l * l / l_0;
  });

  // step 2: add the medial axis energy per node
  func.add_elements<3>(TinyAD::range(active2MedialAxisPairs.size()),
                       [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;
    
    int v0 = active2MedialAxisPairs[e_id].first;
    int v1 = active2MedialAxisPairs[e_id].second;

    Eigen::Vector3<T> x = element.variables(node2ActiveNodeIdx[v0]);

    Eigen::Vector3<T> y = element.variables(node2ActiveNodeIdx[v1]);
    // Eigen::Vector3d y = {nodes[v1].x, nodes[v1].y, nodes[v1].z};

    Eigen::Vector3d medialAxis{active2MedialAxisPositions[e_id].x,
                               active2MedialAxisPositions[e_id].y,
                               active2MedialAxisPositions[e_id].z};

    // === option 1. fix everything except x ===

    // T r = (x - medialAxis).norm();

    // return nodeWeights[v0] * options.alpha * r * r / options.nDirections;

    // === option 2. use the new ln function ===

    Eigen::Vector3d x_ = {nodes[v0].x, nodes[v0].y, nodes[v0].z};
    Eigen::Vector3d y_ = {nodes[v1].x, nodes[v1].y, nodes[v1].z};

    Eigen::Vector3d b = (medialAxis - x_).normalized();
    Eigen::Vector3d xy = y_ - x_;
    double cos_theta = xy.dot(b) / xy.norm();
    double theta = acos(cos_theta);

    T r2 = (x - y).dot(x - y) / cos_theta / 2;
    std::cout << "r2: " << r2 << std::endl;

    return 2*r2 / options.nDirections;
  });

  // func.add_elements<3>(TinyAD::range(active1MedialAxes.size()),
  //                      [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
  //   using T = TINYAD_SCALAR_TYPE(element);
  //   int v_id = active1MedialAxes[element.handle].first;
  //   int activenode_id = node2ActiveNodeIdx[v_id];
  //   auto medialAxis = active1MedialAxes[element.handle].second;

  //   Eigen::Vector3<T> x = element.variables(activenode_id);

  //   Eigen::Vector3d ma {medialAxis.x, medialAxis.y, medialAxis.z};

  //   // === option 1. fix everything ===
  //   T r = (x - ma).norm();

  //   return nodeWeights[v_id] * options.alpha * r * r / options.nDirections;
  // });

  auto x = func.x_from_data([&](int v_idx) {
    auto v = nodes[activeNode2NodeIdx[v_idx]];

    return Eigen::Vector3d(v.x, v.y, v.z);
  });

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  initializing x: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  start = std::chrono::high_resolution_clock::now();

  auto [f, g] = func.eval_with_gradient(x);

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  evaluating f, g, H: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  std::cout << "f: " << f << std::endl;
  // std::cout << "x^TLx: " << x.transpose() * L * x << std::endl;
  Eigen::VectorXd d;

  start = std::chrono::high_resolution_clock::now();

  // project the gradient to the bitangent plane
  // because our energy is invariant under the bitangent plane rotation
  std::vector<Vector3> gradient(nodes.size(), Vector3::zero());

  func.x_to_data(g, [&](int v_idx, const Eigen::Vector3d &p) {
    auto idx = activeNode2NodeIdx[v_idx];

    gradient[idx] = Vector3 {p(0), p(1), p(2)};
  });

  for (int i = 0; i < gradient.size(); i++) {
    if (!isFixedNode[i]) {
      auto v_0 = node2Nodes[i][0];
      auto v_1 = node2Nodes[i][1];

      auto v0 = nodes[v_0];
      auto v1 = nodes[v_1];
      auto v = nodes[i];

      double l0 = (v0 - v).norm();
      double l1 = (v1 - v).norm();

      auto t0 = (v0 - v).normalize();
      auto t1 = (v - v1).normalize();
      auto t = (l0 * t0 + l1 * t1) / (l0 + l1);

      t = t.normalize();

      auto grad = gradient[i];

      gradient[i] = grad - dot(grad, t) * t;
    }
  }

  g = func.x_from_data([&](int v_idx) {
    auto _g = gradient[activeNode2NodeIdx[v_idx]];

    return Eigen::Vector3d(_g.x, _g.y, _g.z);
  });

  Eigen::SimplicialLDLT<SparseMatrix<double>> solver;
  auto [_f, _g, H] = func.eval_with_hessian_proj(x);

  // solver.compute(H);
  // d = -solver.solve(g);
  d = -g;

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  hessian solve: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                    start)
                  .count()
            << " ms" << std::endl;
  
  std::vector<Vector3> descent(nodes.size(), Vector3::zero());

  func.x_to_data(d, [&] (int v_idx, const Eigen::Vector3d& p) {
    auto idx = activeNode2NodeIdx[v_idx];

    descent[idx] = Vector3 {p(0), p(1), p(2)};
  });

  return std::make_tuple(descent, gradient, f);
}
}
