#include "geometrycentral/numerical/linear_algebra_types.h"
#include "modules/volume_filling_energy/newton.h"
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
volume_filling_energy_newton(const std::vector<Vector3> &nodes,
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
        active2MedialAxisPairs.emplace_back(std::pair<int, int> {i, otherNode});
        active2MedialAxisPositions.emplace_back(medialAxisPositions[i][j]);
      } else {
        active1MedialAxes.emplace_back(std::pair<int, Vector3> {i, medialAxisPositions[i][j]});
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

  std::vector<Eigen::Triplet<double>> triplets;
  Eigen::VectorXd b(3*numActiveNodes);

  // step 1: add the curve length energy
  // 2 active nodes
  for (int i = 0; i < segmentsWith2ActiveNodes.size(); i++) {
    int segmentId = activeTwoSegment2SegmentIdx[i];

    int _v0 = segmentsWith2ActiveNodes[i].first;
    int _v1 = segmentsWith2ActiveNodes[i].second;

    int v0 = node2ActiveNodeIdx[_v0];
    int v1 = node2ActiveNodeIdx[_v1];

    double l_0 = segmentLengths[segmentId];

    for (int j = 0; j < 3; j++) {
      int v0id = 3*v0 + j;
      int v1id = 3 * v1 + j;

      triplets.emplace_back(v0id, v0id, options.beta / l_0);
      triplets.emplace_back(v1id, v0id, -options.beta / l_0);
      triplets.emplace_back(v0id, v1id, -options.beta / l_0);
      triplets.emplace_back(v1id, v1id, options.beta / l_0);
    }
  }

  for (int i = 0; i < segmentsWith1ActiveNode.size(); i++) {
    int segmentId = activeOneSegment2SegmentIdx[i];

    int _v0 = segmentsWith1ActiveNode[i].first;
    int _v1 = segmentsWith1ActiveNode[i].second;

    int v0 = node2ActiveNodeIdx[_v0];

    double l_0 = segmentLengths[segmentId];

    for (int j = 0; j < 3; j++) {
      int v0id = 3*v0 + j;

      triplets.emplace_back(v0id, v0id, options.beta / l_0);

      b[v0id] -= nodes[_v1][j] * options.beta / l_0;
    }
  }

  // step 2: add the medial axis energy per node
  for (int i = 0; i < active2MedialAxisPairs.size(); i++) {
    int _v0 = active2MedialAxisPairs[i].first;
    int _v1 = active2MedialAxisPairs[i].second;

    int v0 = node2ActiveNodeIdx[_v0];
    int v1 = node2ActiveNodeIdx[_v1];

    auto medialAxis = active2MedialAxisPositions[i];

    auto x = nodes[_v0];
    auto y = nodes[_v1];
    auto xy = y - x;
    auto bx = (medialAxis - x).normalize();
    auto cos_theta = dot(xy, bx) / xy.norm();
    auto cos_theta_2 = cos_theta * cos_theta;

    for (int j = 0; j < 3; j++) {
      int v0id = 3 * v0 + j;
      int v1id = 3 * v1 + j;

      double coeff = nodeWeights[_v0] * options.alpha / cos_theta_2 / options.nDirections / 8;

      triplets.emplace_back(v0id, v0id, coeff);
      triplets.emplace_back(v0id, v1id, -coeff);
      triplets.emplace_back(v1id, v0id, -coeff);
      triplets.emplace_back(v1id, v1id, coeff);
    }
  }

  for (int i = 0; i < active1MedialAxes.size(); i++) {
    int _v0 = active1MedialAxes[i].first;
    auto medialAxis = active1MedialAxes[i].second;

    int v0 = node2ActiveNodeIdx[_v0];

    for (int j = 0; j < 3; j++) {
      int v0id = 3 * v0 + j;

      triplets.emplace_back(v0id, v0id, nodeWeights[_v0] * options.alpha / options.nDirections);

      b[v0id] -= nodeWeights[_v0] * options.alpha * medialAxis[j] / options.nDirections;
    }
  }

  Eigen::SparseMatrix<double> A(3 * numActiveNodes, 3 * numActiveNodes);
  A.setFromTriplets(triplets.begin(), triplets.end());

  // the energy is approximated as 1/2 * x^T A x + b^T x
  // the hessian is A and the gradient is A x + b

  Eigen::VectorXd x(3 * numActiveNodes);
  for (int i = 0; i < numActiveNodes; i++) {
    auto v = nodes[activeNode2NodeIdx[i]];
    x[3 * i] = v.x;
    x[3 * i + 1] = v.y;
    x[3 * i + 2] = v.z;
  }

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  initializing x: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  start = std::chrono::high_resolution_clock::now();

  double f = (0.5 * x.transpose() * A * x + b.transpose() * x)(0);
  Eigen::VectorXd g = A * x + b;

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  evaluating f, g, H: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  std::cout << "f: " << f << std::endl;
  Eigen::VectorXd d;

  start = std::chrono::high_resolution_clock::now();

  // project the gradient to the bitangent plane
  // because our energy is invariant under the bitangent plane rotation
  std::vector<Vector3> gradient(nodes.size(), Vector3::zero());

  for (int i = 0; i < numActiveNodes; i++) {
    gradient[activeNode2NodeIdx[i]] = Vector3 {g(3 * i), g(3 * i + 1), g(3 * i + 2)};
  }

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

  // auto [_d, _g, _f] = modules::volume_filling_energy_quasi_newton(
  //       nodes, segments, isFixedNode, meshNodes, options);

  for (int i = 0; i < gradient.size(); i++) {
    if (!isFixedNode[i]) {
      int vid = node2ActiveNodeIdx[i];
      g(3 * vid) = gradient[i].x;
      g(3 * vid + 1) = gradient[i].y;
      g(3 * vid + 2) = gradient[i].z;
    }
  }

  start = std::chrono::high_resolution_clock::now();

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

  solver.compute(A);
  d = -solver.solve(g);

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  hessian solve: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                    start)
                  .count()
            << " ms" << std::endl;
  
  std::vector<Vector3> descent(nodes.size(), Vector3::zero());

  for (int i = 0; i < numActiveNodes; i++) {
    descent[activeNode2NodeIdx[i]] =
        Vector3{d(3 * i), d(3 * i + 1), d(3 * i + 2)};
    // std::cout << "descent[activeNode2NodeIdx[i]]: " << descent[activeNode2NodeIdx[i]] << std::endl;
  }

  std::cout << "descent direction computed" << std::endl;

  return std::make_tuple(descent, gradient, f);
}
}
