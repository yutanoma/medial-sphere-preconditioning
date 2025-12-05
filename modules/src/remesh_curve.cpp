#include "modules/remesh_curve.h"

#include "geometrycentral/utilities/vector3.h"
#include <queue>
#include <set>

namespace modules {

// Helper struct to represent an edge with its length
struct Edge {
    int v0, v1;
    double length;
    
    Edge(int v0_, int v1_, double length_) : v0(v0_), v1(v1_), length(length_) {}
    
    bool operator>(const Edge& other) const {
        return length > other.length;
    }
};

std::tuple<
  std::vector<Vector3>,
  std::vector<std::array<int, 2>>,
  std::vector<bool>
> removeShortEdges(
  const std::vector<Vector3> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<bool> &isFixedNode,
  const double &h
) {
    // Priority queue to store edges by length (shortest first)
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> edgeQueue;
    
    // Build adjacency list for quick neighbor lookup
    std::vector<std::vector<int>> adjacencyList(nodes.size());
    std::vector<std::vector<int>> node2Segments(nodes.size());
    
    // Initialize edge queue and adjacency lists
    for (int i = 0; i < segments.size(); i++) {
        int v0 = segments[i][0];
        int v1 = segments[i][1];
        double length = (nodes[v0] - nodes[v1]).norm();
        
        edgeQueue.push(Edge(v0, v1, length));
        
        adjacencyList[v0].push_back(v1);
        adjacencyList[v1].push_back(v0);
        
        node2Segments[v0].push_back(i);
        node2Segments[v1].push_back(i);
    }
    
    // Track deleted nodes and their merging targets
    std::set<int> deletedNodes;
    std::map<int, int> mergeTarget;  // maps deleted node to its target node
    
    // Process edges in order of increasing length
    while (!edgeQueue.empty()) {
        Edge currentEdge = edgeQueue.top();
        edgeQueue.pop();
        
        // Skip if length is above threshold
        if (currentEdge.length >= h) break;
        
        int v0 = currentEdge.v0;
        int v1 = currentEdge.v1;
        
        // Skip if either node has been deleted
        if (deletedNodes.count(v0) > 0 || deletedNodes.count(v1) > 0) continue;
        
        // Skip if both nodes are fixed
        if (isFixedNode[v0] && isFixedNode[v1]) continue;
        
        // Choose which node to keep (prefer fixed nodes)
        int nodeToKeep = isFixedNode[v0] ? v0 : (isFixedNode[v1] ? v1 : v0);
        int nodeToDelete = nodeToKeep == v0 ? v1 : v0;
        
        // Skip if deleting this node would create problematic topology
        if (adjacencyList[nodeToDelete].size() != 2) continue;
        
        // Get the neighbors of the node we're deleting (should be exactly 2)
        std::vector<int> neighbors;
        for (int neighbor : adjacencyList[nodeToDelete]) {
            if (neighbor != nodeToKeep) {  // Skip the edge we're collapsing
                neighbors.push_back(neighbor);
            }
        }
        
        // Mark node as deleted and record merge target
        deletedNodes.insert(nodeToDelete);
        mergeTarget[nodeToDelete] = nodeToKeep;
        
        // Update adjacency list
        adjacencyList[nodeToKeep].clear();  // Clear old adjacencies
        for (int neighbor : neighbors) {
            // Update neighbor's adjacency list
            auto& neighborAdj = adjacencyList[neighbor];
            auto it = std::find(neighborAdj.begin(), neighborAdj.end(), nodeToDelete);
            if (it != neighborAdj.end()) {
                *it = nodeToKeep;
            }
            
            // Add new edge to priority queue
            double newLength = (nodes[nodeToKeep] - nodes[neighbor]).norm();
            edgeQueue.push(Edge(nodeToKeep, neighbor, newLength));
            
            // Update nodeToKeep's adjacency list
            adjacencyList[nodeToKeep].push_back(neighbor);
        }
    }
    
    // Build new node list
    std::vector<Vector3> newNodes;
    std::vector<bool> newNodeIsFixed;
    std::map<int, int> oldToNewIndex;
    
    for (int i = 0; i < nodes.size(); i++) {
        if (deletedNodes.count(i) > 0) continue;
        
        oldToNewIndex[i] = newNodes.size();
        newNodes.push_back(nodes[i]);
        newNodeIsFixed.push_back(isFixedNode[i]);
    }
    
    // Build new segment list
    std::vector<std::array<int, 2>> newSegments;
    for (const auto& segment : segments) {
        int v0 = segment[0];
        int v1 = segment[1];
        
        // Follow merge targets to find final nodes
        while (mergeTarget.count(v0) > 0) v0 = mergeTarget[v0];
        while (mergeTarget.count(v1) > 0) v1 = mergeTarget[v1];
        
        // Skip if both endpoints merged to same node
        if (v0 == v1) continue;
        
        // Add segment with new indices
        newSegments.push_back({oldToNewIndex[v0], oldToNewIndex[v1]});
    }
    
    return {newNodes, newSegments, newNodeIsFixed};
}

std::tuple<
  std::vector<Vector3>,
  std::vector<std::array<int, 2>>,
  std::vector<bool>
> subdivideSegments(
  const std::vector<Vector3> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<bool> &isFixedNode,
  const double &h
) {
  auto newNodes = nodes;
  auto newSegments = segments;
  auto newNodeIsFixed = isFixedNode;

  std::vector<double> segmentLengths(segments.size(), .0);
  for (int i = 0; i < segments.size(); i++) {
    double edgeLen = (nodes[segments[i][0]] - nodes[segments[i][1]]).norm();
    segmentLengths[i] = edgeLen;
  }

  for (int i = 0; i < segments.size(); i++) {
    double edgeLen = segmentLengths[i];

    // if (isFixedNode[segments[i][0]] && isFixedNode[segments[i][1]]) {
    //   continue;
    // }

    if (edgeLen > 2 * h) {
      int divisionNum = std::ceil(edgeLen / (2 * h));
      double lenPerDivision = edgeLen / divisionNum;
      double len = .0;

      std::vector<Vector3> newSurfacePoints = {};

      for (int j = 1; j < divisionNum; j++) {
        len += lenPerDivision;
        double t = len / edgeLen;
        newSurfacePoints.emplace_back((1 - t) * nodes[segments[i][0]] + t * nodes[segments[i][1]]);
      }

      std::vector<int> newNodeIds = {};

      for (int j = 0; j < newSurfacePoints.size(); j++) {
        newNodes.emplace_back(newSurfacePoints[j]);
        newNodeIsFixed.emplace_back(false);
        newNodeIds.emplace_back(newNodes.size() - 1);
      }

      for (int j = 0; j < divisionNum; j++) {
        if (j == 0) {
          newSegments[i][1] = newNodeIds[0];
        } else if (j == divisionNum - 1) {
          newSegments.emplace_back(std::array<int, 2>{newNodeIds[j - 1], segments[i][1]});
        } else {
          newSegments.emplace_back(std::array<int, 2>{newNodeIds[j - 1], newNodeIds[j]});
        }
      }
    }
  }

  return {
    newNodes,
    newSegments,
    newNodeIsFixed
  };
}

std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>, std::vector<bool>>
remesh_curve(const std::vector<Vector3> &nodes,
             const std::vector<std::array<int, 2>> &segments,
             const std::vector<bool> &isFixed,
             double h) {
  // 1. remove nodes that are too close
  auto [newNodes, newSegments, newNodeIsFixed]
    = removeShortEdges(nodes, segments, isFixed, h);

  // 2. subdivide segments
  std::tie(newNodes, newSegments, newNodeIsFixed) = subdivideSegments(newNodes, newSegments, newNodeIsFixed, h);

  if (newNodes.size() < 3 || newSegments.size() < 3) {
    return {
      nodes,
      segments,
      isFixed
    };
  }

  return {
    newNodes,
    newSegments,
    newNodeIsFixed
  };
}
}

