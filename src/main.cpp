#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

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
std::vector<Vector3> meshNodes; // nodes that are fixed for optimization

double h = 0.1;

float timestep = 1.;

bool useGradient = false;
bool run100Iterations = false;
bool run10000Iterations = false;
bool saveCurves = false;

modules::SceneObject scene;
modules::VolumeFillingEnergy::Options options;

int iteration = 0;

std::ofstream timeFile("time.csv");

void doWork() {
  auto prevNodes = nodes;

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

  polyscope::registerCurveNetwork("rest curve", nodes, segments);

  std::vector<Vector3> d;
  std::vector<Vector3> g;
  double f;

  if (scene.preconditioner == "quasi-newton") {
    std::tie(d, g, f) = modules::volume_filling_energy_quasi_newton(
        nodes, segments, isFixed, meshNodes, options);

    polyscope::getPointCloud("points")->addVectorQuantity("gradient", g);
    polyscope::getPointCloud("points")->addVectorQuantity("descent", d, polyscope::VectorType::AMBIENT);
  } else if (scene.preconditioner == "L2") {
    std::tie(d, g, f) = modules::volume_filling_energy_quasi_newton(nodes, segments,
                                                            isFixed, meshNodes, options);
    polyscope::getPointCloud("points")->addVectorQuantity("gradient", g);
    polyscope::getPointCloud("points")->addVectorQuantity(
        "descent", d); //, polyscope::VectorType::AMBIENT);

    
    for (int i = 0; i < nodes.size(); i++) {
      if (!isFixed[i]) {
        d[i] = - h*g[i];
      }
    }
  } else if (scene.preconditioner == "newton") {
    std::tie(d, g, f) = modules::volume_filling_energy_quasi_newton(
        nodes, segments, isFixed, meshNodes, options);

    polyscope::getPointCloud("points")->addVectorQuantity("qn gradient", g);
    polyscope::getPointCloud("points")->addVectorQuantity(
        "qn descent", d); //, polyscope::VectorType::AMBIENT);

    std::tie(d, g, f) = modules::volume_filling_energy_newton(
        nodes, segments, isFixed, meshNodes, options);

    for (int i = 0; i < nodes.size(); i++) {
      if (!isFixed[i]) {
        d[i] = 0.01 * d[i];
      }
    }

    polyscope::getPointCloud("points")->addVectorQuantity(
        "descent", d); //, polyscope::VectorType::AMBIENT);

    polyscope::getPointCloud("points")->addVectorQuantity("gradient", g);


    // auto [_d, _g, _f] = modules::volume_filling_energy_quasi_newton(
    //     nodes, segments, isFixed, meshNodes, options);

    
    // polyscope::getPointCloud("points")->addVectorQuantity(
    //     "descent_quasi_newton", _d); //, polyscope::VectorType::AMBIENT);
    // polyscope::getPointCloud("points")->addVectorQuantity(
    //     "gradient_quasi_newton", _g);

    // std::cout << "f with h1: " << f << std::endl;
    // std::cout << "f with quasi-newton: " << _f << std::endl;

  } else if (scene.preconditioner == "newton_test") {
    std::tie(d, g, f) = modules::volume_filling_energy_newton_test(
        nodes, segments, isFixed, meshNodes, options);

    auto [_d, _g, _f] = modules::volume_filling_energy_quasi_newton(
        nodes, segments, isFixed, meshNodes, options);

    std::cout << "f with h1: " << f << std::endl;
    std::cout << "f with quasi-newton: " << _f << std::endl;

    for (int i = 0; i < nodes.size(); i++) {
      // d[i] = -0.003 * g[i];
      // 
      // d[i] = h * d[i];
    }

    for (int i = 0; i < nodes.size(); i++) {
      // g[i] = h * g[i];
      // _g[i] = h * _g[i];
    }

    polyscope::getPointCloud("points")->addVectorQuantity("gradient", g);
    polyscope::getPointCloud("points")->addVectorQuantity(
        "descent", d); //, polyscope::VectorType::AMBIENT);

    polyscope::getPointCloud("points")->addVectorQuantity(
        "gradient_quasi_newton", _g);
    polyscope::getPointCloud("points")->addVectorQuantity(
        "descent_quasi_newton", _d);

    // for (int i = 0; i < nodes.size(); i++) {
    //   d[i] = 0*d[i];
    // }
  } else if (scene.preconditioner == "laplacian") {
    std::tie(d, g, f) = modules::volume_filling_energy_laplacian(
        nodes, segments, isFixed, meshNodes, options);

    polyscope::getPointCloud("points")->addVectorQuantity("gradient", g);
    polyscope::getPointCloud("points")->addVectorQuantity("descent", d);

  } else {
    std::cerr << "Unknown preconditioner: " << scene.preconditioner << std::endl;
    return;
  }

  std::chrono::steady_clock::time_point ls_begin = std::chrono::steady_clock::now();

  std::tie(nodes, segments, isFixed) = modules::line_search(
      nodes, segments, isFixed, meshNodes, options, d, h, timestep);

  
  std::tie(nodes, segments, isFixed) =
      modules::remesh_curve(nodes, segments, isFixed, h);


  std::chrono::steady_clock::time_point ls_end = std::chrono::steady_clock::now();
  std::cout << "line search done in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(ls_end -
                                                                     ls_begin)
                   .count()
            << " ms" << std::endl;

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "total time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   end - begin)
                   .count()
            << " ms" << std::endl;

  std::cout << "====== iteration " << iteration << " done ======" << std::endl
            << std::endl;

  if (saveCurves) {
    std::string folder = "./curves";
    if (!std::filesystem::exists(folder)) {
      std::filesystem::create_directory(folder);
    }
    std::string filename = folder + "/curve_" + std::to_string(iteration) + ".obj";
    modules::write_curve(filename, nodes, segments);

    timeFile << iteration << "," << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "," << prevNodes.size() << std::endl;
  }

  polyscope::registerCurveNetwork("new curve", nodes, segments);
  iteration++;
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

  ImGui::SliderFloat("alpha", (float*)&options.alpha, 0.0f, 1.0f);
  ImGui::SliderFloat("rMax", (float*)&options.rMax, 0.0f, 1.0f);
  ImGui::SliderInt("nDirections", &options.nDirections, 1, 20);
  ImGui::Checkbox("use gradient", &useGradient);
  ImGui::Checkbox("run 100 iterations", &run100Iterations);
  ImGui::Checkbox("run 10000 iterations", &run10000Iterations);
  ImGui::SliderFloat("timestep", &timestep, 0.0f, 100.0f);
  ImGui::Checkbox("save curves", &saveCurves);

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
  args::Flag excecuteOnly(parser, "executeOnly", "Just excecute the flow for 1000 iterations", {'e'});

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


  timeFile << "iter,time,nodes" << std::endl;

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

  std::cout << "inputFilename: " << args::get(inputFilename) << std::endl;
  std::cout << "rMax: " << scene.rMax << std::endl;
  std::cout << "alpha: " << scene.alpha << std::endl;
  std::cout << "beta: " << scene.beta << std::endl;
  std::cout << "lambda: " << scene.lambda << std::endl;
  std::cout << "h: " << scene.h << std::endl;
  std::cout << "nDirections: " << scene.nDirections << std::endl;

  if (scene.curveFileName == "") {
    std::cerr << "curveFileName: " << scene.curveFileName << " is empty" << std::endl;
    return EXIT_FAILURE;
  }

  std::tie(nodes, segments, isFixed) = modules::read_curve(scene.curveFileName);

  if (scene.meshFileName != "") {
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(scene.meshFileName);

    auto psMesh = polyscope::registerSurfaceMesh(
      "mesh", geometry->inputVertexPositions, mesh->getFaceVertexList());
    psMesh->setTransparency(0.4);

    for (auto v : mesh->vertices()) {
      auto vpos = geometry->inputVertexPositions[v];
      meshNodes.emplace_back(Vector3 {vpos.x, vpos.y, vpos.z});
    }
    polyscope::registerPointCloud("meshNodes", meshNodes)->setEnabled(false);
  }

  std::tie(nodes, segments, isFixed) =
      modules::remesh_curve(nodes, segments, isFixed, h);

  if (excecuteOnly) {
    // create a folder for the curves
    std::string folder = "./curves";
    if (!std::filesystem::exists(folder)) {
      std::filesystem::create_directory(folder);
    }
    for (int i = 0; i < 1000; i++) {
      polyscope::registerPointCloud("points", nodes);
      doWork();
      // save the curves
      std::string filename = "curve_" + std::to_string(i) + ".obj";
      modules::write_curve(folder + "/" + filename, nodes, segments);
    }
    return 0;
  }

  polyscope::registerCurveNetwork("initial curve", nodes, segments)
      ->setEnabled(false);
  polyscope::registerCurveNetwork("rest curve", nodes, segments)->setEnabled(false);
  polyscope::registerPointCloud("points", nodes)->setEnabled(false);

  polyscope::show();

  return 0;
}
