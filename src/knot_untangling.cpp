#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/utilities.h"
#include "imgui.h"
#include "modules/volume_filling_energy/quasi_newton.h"
#include "modules/volume_filling_energy/newton_test.h"
#include "modules/volume_filling_energy/newton.h"
#include "modules/volume_filling_energy/laplacian.h"
#include "modules/line_search.h"
#include "modules/remesh_curve.h"
#include "modules/scene_file.h"
#include "modules/read_curve.h"
#include "modules/write_curve.h"
#include "modules/medial_axis_3d.h"
#include "polyscope/affine_remapper.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <iostream>

#include <igl/triangle/triangulate.h>

#include "args/args.hxx"

using namespace geometrycentral;
using namespace geometrycentral::surface;

std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

std::vector<Vector3> nodes;
std::vector<std::array<int, 2>> segments;
std::vector<bool> isFixed;

double h = 0.1;

double original_h = 0;
double original_alpha = 0;

float timestep = 1.;
double originalLength = 0;

bool useGradient = false;
bool run100Iterations = false;
bool run10000Iterations = false;
bool saveCurves = false;
bool normalizeLength = false;
modules::SceneObject scene;
modules::VolumeFillingEnergy::Options options;

int iteration = 0;

std::ofstream gradNormFile("gradNorm.txt");
std::ofstream timeFile("time.txt");

void doWork() {
  std::cout << "====== iteration " << iteration << " start ======" << std::endl;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::cout << "nodes: " << nodes.size() << std::endl;
  std::cout << "segments: " << segments.size() << std::endl;
  int fixedCount = 0;
  for (int i = 0; i < isFixed.size(); i++) {
    if (isFixed[i]) {
      fixedCount++;
    }
  }
  std::cout << "fixed nodes: " << fixedCount << std::endl;

  double oldLength = 0;
  for (int i = 0; i < segments.size(); i++) {
    auto v0 = nodes[segments[i][0]];
    auto v1 = nodes[segments[i][1]];
    oldLength += (v0 - v1).norm();
  }
  std::cout << "old length: " << oldLength << std::endl;

  // h = oldLength / 600;

  std::vector<Vector3> d;
  std::vector<Vector3> g;
  double f;

  std::vector<Vector3> meshNodes;
  std::chrono::steady_clock::time_point begin_preconditioner = std::chrono::steady_clock::now();
  if (scene.preconditioner == "laplacian") {
    std::tie(d, g, f) = modules::volume_filling_energy_laplacian(
        nodes, segments, isFixed, meshNodes, options);
  } else if (scene.preconditioner == "quasi-newton") {
    std::tie(d, g, f) = modules::volume_filling_energy_quasi_newton(
        nodes, segments, isFixed, meshNodes, options);
  } else if (scene.preconditioner == "newton") {
    std::tie(d, g, f) = modules::volume_filling_energy_newton(
        nodes, segments, isFixed, meshNodes, options);
  } else {
    throw std::runtime_error("preconditioner not supported");
  }
  std::chrono::steady_clock::time_point end_preconditioner = std::chrono::steady_clock::now();
  std::cout << "preconditioner time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   end_preconditioner - begin_preconditioner)
                   .count()
            << " ms" << std::endl;
  // for (int i = 0; i < nodes.size(); i++) {
  //   d[i] = h*d[i];
  // }

  polyscope::registerCurveNetwork("rest curve", nodes, segments);

  std::chrono::steady_clock::time_point ls_begin = std::chrono::steady_clock::now();

  std::tie(nodes, segments, isFixed) =
      modules::line_search(nodes, segments, isFixed, meshNodes, options, d, h, timestep);
      

  std::chrono::steady_clock::time_point ls_end = std::chrono::steady_clock::now();
  std::cout << "line search done in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(ls_end -
                                                                     ls_begin)
                   .count()
            << " ms" << std::endl;

    double newLength = 0;
  for (int i = 0; i < segments.size(); i++) {
    auto v0 = nodes[segments[i][0]];
    auto v1 = nodes[segments[i][1]];
    newLength += (v0 - v1).norm();
  }
  std::cout << "new length: " << newLength << std::endl;

  double gradNorm = 0;
  for (int i = 0; i < g.size(); i++) {
    gradNorm += g[i].norm();
  }
  double r = 3.36 / std::sqrt(options.alpha);
  double gradNormRatio = gradNorm / nodes.size() / r;
  std::cout << "gradNorm: " << gradNormRatio << std::endl;
  // gradNormFile << gradNorm / nodes.size() << std::endl;
  gradNormFile << gradNormRatio << std::endl;



  if (normalizeLength) {
    std::vector<Vector3> scaledNodes;
    Vector3 gc = Vector3 {0, 0, 0};
    for (int i = 0; i < nodes.size(); i++) {
      scaledNodes.emplace_back(nodes[i] * originalLength / newLength);
      gc += scaledNodes[i];
    }
    gc /= nodes.size();
    for (int i = 0; i < nodes.size(); i++) {
      scaledNodes[i] -= gc;
    }
    nodes = scaledNodes;
    double r = originalLength / M_PI / 2;
    // options.alpha = 3.36 * 3.36 / r / r;
    // options.alpha = 0;
    // options.rMax = 0.00001;
    // options.rMax = r - 1 / r / 2;
    // options.rMax = 0.9 * r;
    // options.rMax = 0.00001;
    // h = r / 5;
    std::cout << "r: " << r << std::endl;
    std::cout << "alpha: " << options.alpha << std::endl;
    std::cout << "rMax: " << options.rMax << std::endl;
    std::cout << "h: " << h << std::endl;
  }

  std::tie(nodes, segments, isFixed) =
      modules::remesh_curve(nodes, segments, isFixed, h);

  std::vector<Vector3> scaledNodes;
  Vector3 gc = Vector3 {0, 0, 0};
  for (int i = 0; i < nodes.size(); i++) {
    scaledNodes.emplace_back(nodes[i] * originalLength / newLength);
    gc += scaledNodes[i];
  }
  gc /= nodes.size();
  for (int i = 0; i < nodes.size(); i++) {
    scaledNodes[i] -= gc;
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "total time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   end - begin)
                   .count()
            << " ms" << std::endl;
  timeFile << std::chrono::duration_cast<std::chrono::milliseconds>(
                   end - begin)
                   .count() << std::endl;

  std::cout << "====== iteration " << iteration << " done ======" << std::endl
            << std::endl;

  polyscope::getPointCloud("points")->addVectorQuantity("gradient", g);
  polyscope::getPointCloud("points")->addVectorQuantity("descent", d);

  polyscope::registerCurveNetwork("new curve", scaledNodes, segments);
  iteration++;

  if (saveCurves) {
    std::string folder = "./curves";
    if (!std::filesystem::exists(folder)) {
      std::filesystem::create_directory(folder);
    }
    std::string filename = folder + "/curve_" + std::to_string(iteration) + ".obj";
    modules::write_curve(filename, nodes, segments);

    auto edgetails = polyscope::getCurveNetwork("maPairs")->edgeTailInds;
    auto edgetips = polyscope::getCurveNetwork("maPairs")->edgeTipInds;
    auto nodepositions = polyscope::getCurveNetwork("maPairs")->nodePositions;
    std::ofstream ofs("./curves/maPairs_" + std::to_string(iteration) + ".obj");
    for (int i = 0; i < edgetails.size(); i++) {
      auto v0 = nodepositions.getValue(edgetails.getValue(i));
      auto v1 = nodepositions.getValue(edgetips.getValue(i));
      ofs << "v " << v0.x << " " << v0.y << " " << v0.z << std::endl;
      ofs << "v " << v1.x << " " << v1.y << " " << v1.z << std::endl;
    }
    for (int i = 0; i < edgetails.size(); i++) {
      ofs << "l " << 2 * i + 1 << " " << 2 * i + 2 << std::endl;
    }
  }

  if (gradNormRatio < 0.005 || nodes.size() < 7) {
    std::string filename = "./curves/curve_final.obj";
    modules::write_curve(filename, nodes, segments);
    std::cout << "converged" << std::endl;
    std::exit(0);
  }
}

int ranIterations = 0;

void myCallback() {
  if (ImGui::IsKeyDown(' ') || run100Iterations || run10000Iterations) {
    polyscope::registerPointCloud("points", nodes);
    doWork();
  }

  if (ImGui::Button("compute one step")) {
    polyscope::registerPointCloud("points", nodes);
    doWork();
  }

  ImGui::Checkbox("run 100 iterations", &run100Iterations);
  ImGui::Checkbox("run 10000 iterations", &run10000Iterations);
  ImGui::Checkbox("save curves", &saveCurves);
  ImGui::Checkbox("normalize length", &normalizeLength);
  ImGui::SliderFloat("timestep", &timestep, 0.0f, 100.0f);
  if (run100Iterations || run10000Iterations) {
    ranIterations++;
    if (run100Iterations && ranIterations == 100) {
      run100Iterations = false;
      ranIterations = 0;
    }
    if (run10000Iterations && ranIterations == 10000) {
      run10000Iterations = false;
      ranIterations = 0;
    }
  }
}

int main(int argc, char **argv) {
  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
  args::Flag excecuteOnly(parser, "excecuteOnly", "Excecute only", {"excecuteOnly", "e"});
  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a scene file was given
  if (!inputFilename) {
    std::cerr << "Please specify a scene file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Delete all files in ./curves if exists, and create the folder if not
  std::string folder = "./curves";
  if (std::filesystem::exists(folder)) {
    for (const auto& entry : std::filesystem::directory_iterator(folder)) {
      std::filesystem::remove_all(entry.path());
    }
  } else {
    std::filesystem::create_directory(folder);
  }


  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load scene
  scene = modules::read_scene(args::get(inputFilename));


  options.rMax = scene.rMax;
  options.alpha = scene.alpha;
  options.beta = scene.beta;
  options.lambda = scene.lambda;
  options.nDirections = scene.nDirections;
  h = scene.h;
  timestep = scene.timestep;

  double r = 3.36 / std::sqrt(scene.alpha);

  std::cout << "inputFilename: " << args::get(inputFilename) << std::endl;
  std::cout << "rMax: " << scene.rMax << std::endl;
  std::cout << "alpha: " << scene.alpha << std::endl;
  std::cout << "h: " << scene.h << std::endl;
  std::cout << "estimated r: " << r << std::endl;

  if (scene.curveFileName == "") {
    std::cerr << "curveFileName: " << scene.curveFileName << " is empty" << std::endl;
    return EXIT_FAILURE;
  }

  std::tie(nodes, segments, isFixed) = modules::read_curve(scene.curveFileName);

  std::tie(nodes, segments, isFixed) =
      modules::remesh_curve(nodes, segments, isFixed, h);

  polyscope::registerPointCloud("points", nodes);

  std::vector<Vector3> meshNodes;
  auto [medialAxisPositions, otherNodes] = modules::medial_axis_3d(
      nodes, segments, meshNodes, 1, options.nDirections);

  if (!scene.useVal) {
    double r_min = 1e10;
    double r_avg = 0;
    int otherNodesCount = 0;
    for (int i = 0; i < medialAxisPositions.size(); i++) {
      for (int j = 0; j < medialAxisPositions[i].size(); j++) {
        r_min = std::min(r_min, (nodes[i] - medialAxisPositions[i][j]).norm());
        if (otherNodes[i][j] != -1) {
          r_avg += (nodes[i] - medialAxisPositions[i][j]).norm();
          otherNodesCount++;
        }
      }
    }
    r_avg /= otherNodesCount;
    std::cout << "r_min: " << r_min << std::endl;
    std::cout << "r_avg: " << r_avg << std::endl;

    options.rMax = r_avg * 0.9;
    options.alpha = std::pow(3.36, 2) / std::pow(r_avg, 2);
    r = 3.36 / std::sqrt(options.alpha);
    h = r / 5;
    std::cout << std::endl;
    std::cout << "rMax: " << options.rMax << std::endl;
    std::cout << "alpha: " << options.alpha << std::endl;
    std::cout << "h: " << h << std::endl;
    std::cout << "estimated r: " << r << std::endl; 
  }

  double oldLength = 0;
  for (int i = 0; i < segments.size(); i++) {
    auto v0 = nodes[segments[i][0]];
    auto v1 = nodes[segments[i][1]];
    oldLength += (v0 - v1).norm();
  }
  originalLength = oldLength;

  original_h = h;
  original_alpha = options.alpha;

  std::tie(nodes, segments, isFixed) = modules::remesh_curve(nodes, segments, isFixed, h);

  polyscope::registerCurveNetwork("initial curve", nodes, segments)
      ->setEnabled(false);
  polyscope::registerCurveNetwork("rest curve", nodes, segments)->setEnabled(false);
  polyscope::registerPointCloud("points", nodes)->setEnabled(false);

  polyscope::registerCurveNetwork("new curve", nodes, segments);

  if (excecuteOnly) {
    std::string folder = "./curves";
    if (!std::filesystem::exists(folder)) {
      std::filesystem::create_directory(folder);
    }
    for (int i = 0; i < 10000; i++) {
      polyscope::registerPointCloud("points", nodes);
      doWork();
      std::string filename = "./curves/curve_" + std::to_string(i) + ".obj";
      modules::write_curve(filename, nodes, segments);

      if (i == 9999) {
        std::cout << "not converged" << std::endl;
        std::exit(1);
      }
    }
  }

  polyscope::show();

  return 0;
}
