#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <igl/PI.h>

namespace modules {
  enum CurveType {
    Closed,
    Open
  };

  struct SceneObject {
    std::string curveFileName = "";
    std::string meshFileName = "";
    double rMax = 1;
    double alpha = 10;
    double beta = 0.01;
    double lambda = 1;
    double h = 0.1;
    int nDirections = 32;
    std::string preconditioner = "laplacian";
    double timestep = 1;
    bool useVal = false;
  };

  void splitString(const std::string &str, std::vector<std::string> &cont, char delim = ' ');

  SceneObject read_scene(std::string filename);
}
