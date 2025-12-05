#include "modules/read_curve.h"

#include "modules/scene_file.h"

const double EPSILON = 1e-5;

namespace modules {
std::tuple < std::vector<Vector3>,
    std::vector<std::array<int, 2>>,
    std::vector<bool>
>
read_curve(const std::string filename) {
  std::ifstream inFile;
  inFile.open(filename);

  if (!inFile) {
    std::cerr << "Could not open file " << filename << std::endl;
    exit(1);
  }

  std::vector<Vector3> addedNodes = {};
  std::map<int, bool> isAddedNodeFixed = {};
  std::vector<std::array<int, 2>> addedSegments = {};

  for (std::string line; std::getline(inFile, line ); ) {
    std::vector<std::string> parts;
    modules::splitString(line, parts, ' ');

    if (parts.size() == 0) continue;

    if (parts[0] == "v" && parts.size() >= 4) {
      auto coord = Vector3 {
        std::stof(parts[1]),
        std::stof(parts[2]),
        std::stof(parts[3])
      };

      int vid = addedNodes.size();
      addedNodes.push_back(coord);
    } else if (parts[0] == "l") {
      for (int i = 2; i < parts.size(); i++) {
        int v0 = (std::stoi(parts[i-1]) - 1);
        int v1 = (std::stoi(parts[i]) - 1);
        std::array<int, 2> line = {v0, v1};
        addedSegments.push_back(line);
      }
    }
    else if (parts[0] == "fixed") {
      for (int i = 1; i < parts.size(); i++) {
        int id = (std::stoi(parts[i]) - 1);
        isAddedNodeFixed[id] = true;
      }
    }
  }

  std::vector<bool> isFixed(addedNodes.size(), false);
  for (const auto& [id, fixed] : isAddedNodeFixed) {
    isFixed[id] = fixed;
    std::cout << "fixed node " << id << " is " << fixed << std::endl;
  }

  return {addedNodes, addedSegments, isFixed};
}
}
