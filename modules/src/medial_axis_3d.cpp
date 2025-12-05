#include "modules/medial_axis_3d.h"

#include <modules/knn-cpp/knncpp.h>

#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

namespace modules {
int closestPointIndex(
  const Vector3 &p,
  const knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> &kdtree
) {
  Eigen::MatrixXd P(3, 1);
  Eigen::MatrixXd distances;
  knncpp::Matrixi indices;
  P << p.x, p.y, p.z;
  kdtree.query(P, 1, indices, distances);

  assert(indices.rows() == 1);

  return indices(0, 0);
}

std::tuple<double, int> maximumBallRadius(
  const Vector3 &x,
  const Vector3 &b,
  const int i,
  const double maxRadius,
  const Eigen::MatrixXd &cartesianCoords,
  knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> &kdtree
) {
  auto r = maxRadius;
  auto c = x + r * b;

  int nn = closestPointIndex(c, kdtree);
  bool finished = nn == i;

  double bsMax = 1., bsMin = 0.;
  int itrc = 0;

  // continue until the radius is found
  while (!finished) {
    // std::cout << "  " << bsMax << ", " << bsMin << ", " << (bsMax + bsMin) / 2 << std::endl;

    itrc++;

    r = maxRadius * (bsMax + bsMin) / 2;

    auto c = x + r * b;
    int nn = closestPointIndex(c, kdtree);

    if (nn == i) {
      bsMin = (bsMax + bsMin) / 2;
    } else {
      // check if there is a point inside the ball
      auto nn_vec = cartesianCoords.col(nn);
      auto i_vec = cartesianCoords.col(i);
      auto xy_vec = nn_vec - i_vec;
      Vector3 xy = Vector3{xy_vec(0), xy_vec(1), xy_vec(2)};
      r = pow(xy.norm(), 2) / (2 * dot(xy, b));

      auto c = x + r * b;
      int _nn = closestPointIndex(c, kdtree);

      auto _nn_vec = cartesianCoords.col(_nn);
      Vector3 _nn_pos = Vector3{_nn_vec(0), _nn_vec(1), _nn_vec(2)};

      double dist = (_nn_pos - c).norm();
      double error = std::abs(dist - r);

      if (_nn == nn || _nn == i || error < 1e-6) {
        // if no one inside, break
        finished = true;
        return {r, nn};
      } else {
        // else, continue the binary search
        bsMax = (bsMax + bsMin) / 2;
        assert(bsMax > bsMin);
      }
    }

    if (itrc > 100) {
      break;
    }
  }

  // std::cout << itrc << std::endl;

  // could not find a radius
  return {r, -1};
}

std::vector<Vector3>
tangent_vectors(const std::vector<Vector3> &nodes,
                const std::vector<std::array<int, 2>> &segments) {
  std::vector<Vector3> tangents(nodes.size());
  
  // Construct vertex to segment adjacency
  std::vector<std::vector<size_t>> vertexToSegments(nodes.size());
  for (size_t i = 0; i < segments.size(); i++) {
    vertexToSegments[segments[i][0]].push_back(i);
    vertexToSegments[segments[i][1]].push_back(i);
  }

  // Compute tangents for each vertex
  for (size_t v = 0; v < nodes.size(); v++) {
    const auto& adjSegments = vertexToSegments[v];
    
    if (adjSegments.empty()) {
      // No segments connected to this vertex
      tangents[v] = Vector3::zero();
      continue;
    }
    
    if (adjSegments.size() == 1) {
      // Only one segment, use its direction
      size_t segIdx = adjSegments[0];
      int v0 = segments[segIdx][0];
      int v1 = segments[segIdx][1];
      Vector3 t = (nodes[v1] - nodes[v0]).normalize();
      // Flip direction if this vertex is the end point
      if (v == v1) t = -t;
      tangents[v] = t;
      continue;
    }
    
    // For vertices with multiple segments, use first two segments
    // and ensure proper directionality
    Vector3 avgTangent = Vector3::zero();
    for (size_t i = 0; i < std::min(size_t(2), adjSegments.size()); i++) {
      size_t segIdx = adjSegments[i];
      int v0 = segments[segIdx][0];
      int v1 = segments[segIdx][1];

      // weight the tangent vector by the length of the segment
      double length = (nodes[v1] - nodes[v0]).norm();

      // let's say the segment looks like u0 -- v -- u1
      // we want to compute the tangent vector for v
      // the first segment is u0 -> v and the second is v -> u1
      int u = v0 == v ? v1 : v0;
      if (i == 0) {
        // if i == 0 we want to get the direction of u -> v
        avgTangent += (nodes[v] - nodes[u]).normalize() * length;
      } else if (i == 1) {
        // if i == 1 we want to get the direction of v -> u
        avgTangent += (nodes[u] - nodes[v]).normalize() * length;
      }
    }
    tangents[v] = avgTangent.normalize();
  }

  return tangents;
}

std::tuple<std::vector<std::vector<Vector3>>, std::vector<std::vector<int>>>
medial_axis_3d(const std::vector<Vector3> &nodes,
               const std::vector<std::array<int, 2>> &segments,
               const std::vector<Vector3> &meshNodes,
               const double maxRadius,
               const int nDirections) {
    // create a closest point query
    Eigen::MatrixXd _V(3, nodes.size() + meshNodes.size());
    for (int i = 0; i < nodes.size(); i++) {
      _V.col(i) << nodes[i].x, nodes[i].y, nodes[i].z;
    }

    // add the mesh nodes
    // this is only added to the closest point query and not used for anything else
    for (int i = 0; i < meshNodes.size(); i++) {
      _V.col(nodes.size() + i) << meshNodes[i].x, meshNodes[i].y, meshNodes[i].z;
    }

    knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> kdtree(_V);
    kdtree.build();

    std::vector<std::vector<Vector3>> nodeMedialAxis(nodes.size());
    std::vector<std::vector<int>> otherNode(nodes.size());

    auto nodeTangents = tangent_vectors(nodes, segments);
    std::vector<std::vector<Vector3>> rayDirections(nDirections);

    for (int i = 0; i < nodes.size(); i++) {
      auto t = nodeTangents[i].normalize();
      auto x = nodes[i];

      Vector3 b{1, 0, 0};

      if (std::abs(dot(t, b)) > 0.999) {
        b = Vector3{0, 1, 0};
      }

      Vector3 n = cross(t, b).normalize();
      b = cross(n, t).normalize();

      for (int j = 0; j < nDirections; j++) {
        // rotate b around t by 2 * pi / nDirections
        // Construct rotation matrix around axis t by angle theta using Rodrigues' formula
        double theta = 2 * M_PI / nDirections;
        Eigen::Vector3d axis(t.x, t.y, t.z);
        Eigen::Vector3d v(b.x, b.y, b.z);
        
        // Rodrigues' rotation formula
        Eigen::Matrix3d K;  // cross-product matrix
        K << 0, -axis(2), axis(1),
             axis(2), 0, -axis(0),
             -axis(1), axis(0), 0;
             
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + 
                           std::sin(theta) * K +
                           (1 - std::cos(theta)) * K * K;
                           
        Eigen::Vector3d rotated = R * v;
        b = Vector3{rotated(0), rotated(1), rotated(2)};

        // std::cout << "[" << i << ", " << j << "]: " << b << std::endl;

        auto [r_min, nn] = maximumBallRadius(x, b, i, maxRadius, _V, kdtree);
        if (nn == -1 || r_min > maxRadius) {
          r_min = maxRadius;
          nn = -1;
        } else if (nn >= nodes.size()) {
          nn = -1;
        }

        // std::cout << "r_min: " << r_min << " " << maxRadius << std::endl;
        // std::cout << "b: " << b << std::endl;

        // std::cout << "result: " << r_min_minus << ", " << r_min_plus <<
        // std::endl;
        // std::cout << std::endl;

        nodeMedialAxis[i].emplace_back(nodes[i] + r_min * b);
        otherNode[i].emplace_back(nn);
        rayDirections[j].emplace_back(b);
      }
    }

    for (int i = 0; i < nDirections; i++) {
      polyscope::getPointCloud("points")->addVectorQuantity("rayDirections" + std::to_string(i), rayDirections[i]);
    }

    return {nodeMedialAxis, otherNode};
}
} // namespace modules
