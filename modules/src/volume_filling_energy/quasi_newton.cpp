#include "modules/volume_filling_energy/quasi_newton.h"
#include "modules/medial_axis_3d.h"

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <cstddef>
#include <fstream>
#include <map>

#include <chrono>

#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

namespace modules {
std::tuple<std::vector<Vector3>, std::vector<Vector3>, double>
volume_filling_energy_quasi_newton(
    const std::vector<Vector3> &nodes,
    const std::vector<std::array<int, 2>> &segments,
    const std::vector<bool> &isFixedNode, const std::vector<Vector3> &meshNodes,
    const VolumeFillingEnergy::Options &options,
    const bool energyOnly) {
  std::cout << "ndirections: " << options.nDirections << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<int>> node2Segments(nodes.size(), std::vector<int>{});
  for (int i = 0; i < segments.size(); i++) {
    node2Segments[segments[i][0]].emplace_back(i);
    node2Segments[segments[i][1]].emplace_back(i);
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

  // for (int i = 0; i < medialAxisPositions.size(); i++) {
  //   for (int j = 0; j < medialAxisPositions[i].size(); j++) {
  //     if (otherNodes[i][j] == -1) {
  //       auto node = nodes[i];
  //       auto medialAxis = medialAxisPositions[i][j];
  //       auto b = (node - medialAxis).normalize();
  //       medialAxisPositions[i][j] = medialAxis + options.rMax * b / 4;
  //     }
  //   }
  // }

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  sampling medial axis on S^n: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  start = std::chrono::high_resolution_clock::now();

  {
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

    polyscope::registerCurveNetwork("maPairs", nodes, maPairs);

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

  std::vector<double> alphas(numActiveNodes, options.alpha);

  if (options.varyingAlpha) {
    if (options.varyingAlphaType == "y") {
      // 1. gradient according to z
      std::vector<double> yAxis(numActiveNodes, 0.0);
      double miny = 1e10, maxy = -1e10;
      for (int j = 0; j < numActiveNodes; j++) {
        yAxis[j] = nodes[activeNode2NodeIdx[j]].y;
        miny = std::min(miny, yAxis[j]);
        maxy = std::max(maxy, yAxis[j]);
      }

      double maxAlpha = 4 * options.alpha;
      double minAlpha = options.alpha;

      for (int j = 0; j < numActiveNodes; j++) {
        alphas[j] = minAlpha + (maxAlpha - minAlpha) * (yAxis[j] - miny) / (maxy - miny);
      }
    } else if (options.varyingAlphaType == "circle") {
      // 1. get the barycenter of the active nodes
      Vector3 barycenter = Vector3::zero();
      for (int j = 0; j < numActiveNodes; j++) {
        barycenter += nodes[activeNode2NodeIdx[j]];
      }
      barycenter /= numActiveNodes;

      // 2. get the radius of the circle
      double maxRadius = 0.0;
      for (int j = 0; j < numActiveNodes; j++) {
        maxRadius = std::max(maxRadius, (nodes[activeNode2NodeIdx[j]] - barycenter).norm());
      }

      // 3. get the alpha for each node
      for (int j = 0; j < numActiveNodes; j++) {
        double distance = (nodes[activeNode2NodeIdx[j]] - barycenter).norm();
        if (distance > 0.1) {
          alphas[j] = options.alpha;
        } else {
          alphas[j] = options.alpha * 4;
        }
      }
    } else if (options.varyingAlphaType == "gaussian") {
      
    }
  }

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

    return l * l / l_0;
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

    return l * l / l_0;
  });

  // step 2: add the medial axis energy per node
  // in quasi-newton, we can fix the medial axis positions and treat it as a constant
  func.add_elements<1>(TinyAD::range(numActiveNodes),
    [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
      using T = TINYAD_SCALAR_TYPE(element);
      int v_id = element.handle;

      Eigen::Vector3<T> x = element.variables(v_id);
      auto medialAxis =
          medialAxisPositions[activeNode2NodeIdx[v_id]];

      T medialAxisEnergy = 0.0;

      for (int i = 0; i < medialAxis.size(); i++) {
        auto _y = medialAxis[i];
        auto y = Eigen::Vector3d{_y.x, _y.y, _y.z};
        auto r = (x - y).norm();

        // if (otherNodes[activeNode2NodeIdx[v_id]][i] == -1) {
        //   medialAxisEnergy += r * r * 3.6;
        // } else {
        //   medialAxisEnergy += r * r;
        // }
        medialAxisEnergy += r * r;
      }

      // std::cout << "qn medialAxisEnergy [" <<
      // activeNode2NodeIdx[v_id]
      //           << "]: " << medialAxisEnergy << std::endl;
      // std::cout << "activenode: " << v_id << std::endl;
      // std::cout << "node:"
      // std::cout << "options.nDirections: " <<
      // options.nDirections << std::endl;

      double alpha = alphas[activeNode2NodeIdx[v_id]];

      return nodeWeights[activeNode2NodeIdx[v_id]] *
            alpha * medialAxisEnergy /
            options.nDirections;
    });

  auto x = func.x_from_data([&](int v_idx) {
    auto v = nodes[activeNode2NodeIdx[v_idx]];

    return Eigen::Vector3d(v.x, v.y, v.z);
  });

  // for (int i = 0; i < numActiveNodes; i++) {
  //   std::cout << "input node: v[" << i << "]: " << nodes[activeNode2NodeIdx[i]] << std::endl;
  // }
  // for (int i = 0; i < x.rows(); i++) {
  //   std::cout << "output node: x[" << i << "]: " << x[i] << std::endl;
  // }

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  initializing x: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  start = std::chrono::high_resolution_clock::now();

  if (energyOnly) {
    double f = func.eval(x);
    return std::make_tuple(std::vector<Vector3>(), std::vector<Vector3>(), f);
  }

  TinyAD::LinearSolver solver;
  auto [f, g, H] = func.eval_with_hessian_proj(x);

  end = std::chrono::high_resolution_clock::now();
  std::cout << "  evaluating f, g, H: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << " ms" << std::endl;
  Eigen::VectorXd d;


  start = std::chrono::high_resolution_clock::now();

  d = TinyAD::newton_direction(g, H, solver);

  auto x_new = TinyAD::line_search(x, d, f, g, func);
  d = x_new - x;


  end = std::chrono::high_resolution_clock::now();
  std::cout << "  hessian solve: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                    start)
                  .count()
            << " ms" << std::endl;
  
  std::vector<Vector3> descent(nodes.size(), Vector3::zero());
  std::vector<Vector3> gradient(nodes.size(), Vector3::zero());

  func.x_to_data(d, [&] (int v_idx, const Eigen::Vector3d& p) {
    auto idx = activeNode2NodeIdx[v_idx];

    descent[idx] = Vector3 {p(0), p(1), p(2)};
  });

  func.x_to_data(g, [&](int v_idx, const Eigen::Vector3d &p) {
    auto idx = activeNode2NodeIdx[v_idx];

    gradient[idx] = Vector3 {p(0), p(1), p(2)};
  });

  {
    // std::ofstream ofs("./qn-preconditioner.txt");
    // for (int i = 0; i < H.rows(); i++) {
    //   for (int j = 0; j < H.cols(); j++) {
    //     if (H.coeff(i, j) != 0) {
    //       ofs << i << " " << j << " " << H.coeff(i, j) << std::endl;
    //     }
    //   }
    // }
  }

  // {
  //   std::ofstream file("radius.txt");
  //   for (int i = 0; i < nodes.size(); i++) {
  //     for (int j = 0; j < medialAxisPositions[i].size(); j++) {
  //       auto r = (nodes[i] - medialAxisPositions[i][j]).norm();
  //       file << r << std::endl;
  //     }
  //   }
  //   file.close();
  // }

  return std::make_tuple(descent, gradient, f);
}
}
