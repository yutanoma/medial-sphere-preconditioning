#include "modules/scene_file.h"

namespace modules {
  using namespace std;

  void splitString(const std::string& str, std::vector<string> &cont, char delim) {
      std::stringstream ss(str);
      std::string token;
      while (std::getline(ss, token, delim)) {
        cont.push_back(token);
      }
  }

  std::string getDirectoryFromPath(std::string str) {
    using namespace std;
    vector<string> parts;
    splitString(str, parts, '/');

    int nParts = parts.size();
    if (nParts == 1) return "./";
    
    string path = "";

    for (int i = 0; i < nParts - 1; i++) {
        path = path + parts[i] + "/";
    }

    return path;
  }

  void processLine(SceneObject &scene, std::string directory, std::vector<std::string> &parts) {
    string key = parts[0];

    if (key == "#" || key == "//") {
      return;
    }

    if (key == "curve") {
      scene.curveFileName = directory + parts[1];
    } else if (key == "mesh") {
      scene.meshFileName = directory + parts[1];
    } else if (key == "curve") {
      scene.meshFileName = directory + parts[1];
    } else if (key == "rmax") {
      scene.rMax = stod(parts[1]);
    } else if (key == "alpha") {
      scene.alpha = stod(parts[1]);
    } else if (key == "beta") {
      scene.beta = stod(parts[1]);
    } else if (key == "lambda") {
      scene.lambda = stod(parts[1]);
    } else if (key == "h") {
      scene.h = stod(parts[1]);
      std::cout << "h: " << scene.h << std::endl;
    } else if (key == "nDirections") {
      scene.nDirections = stoi(parts[1]);
    } else if (key == "preconditioner") {
      scene.preconditioner = parts[1];
    } else if (key == "timestep") {
      scene.timestep = stod(parts[1]);
    } else if (key == "useVal") {
      scene.useVal = true;
    } else {
      std::cerr << "Unknown key: " << key << std::endl;
    }
  }

  SceneObject read_scene(std::string filename) {
    string directory = getDirectoryFromPath(filename);

    ifstream inFile;
    inFile.open(filename);

    if (!inFile) {
      cerr << "Could not open file " << filename << endl;
      exit(1);
    }

    SceneObject scene;

    std::vector<std::string> parts;
    for (std::string line; std::getline(inFile, line ); ) {
      if (line == "" || line == "\n") continue;
      parts.clear();
      splitString(line, parts, ' ');
      processLine(scene, directory, parts);
    }

    inFile.close();
    return scene;
  }
}
