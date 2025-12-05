## Medial Axis Preconditioning for Knot Untangling and Volume-Filling Curves [SIGGRAPH Asia 2025]

![Medial Axis Preconditioning Teaser](teaser.jpg)

Public code release for Medial Axis Preconditioning for Knot Untangling and Volume-Filling Curves.

```
@article{Noma2025MedialSphere,
  title = {Medial Sphere Preconditioning for Knot Untangling and Volume-Filling Curves},
  author = {Yuta Noma and Alec Jacobson and Karan Singh},
  year = {2025},
  journal = {SIGGRAPH Asia Conference Papers}, 
}
```

## Build the code

We require that `boost` is installed in the user's system.

First, install all the dependencies with `git submodule update --init --recursive`.

### Unix-like machines
Configure (with cmake) and compile.

```
cd /path/to/directory
mkdir build
cd build
cmake ..
make -j
```

### Windows / Visual Studio

Install CMake, and use either the CMake GUI or the command line interface (as on unix) to generate a Visual Studio solution. Build the solution with Visual Studio.

## Run the code

Building will produce two binaries: `main` and `knot_untanling`. `main` produces a volume-filling curve using the given parameters, while `knot_untangling` untangles a given knot using the hyperparameters we provided in Section 4.4 in the paper.

For either case, please press the space bar to run the geometric flow.

### Volume-filling curves

```
./bin/main ../models/bunny/scene.txt
```

### Knot untangling

```
./bin/knot_untangling ../models/ochiai1/scene.txt
```

### Contact

If any issues or questions, please contact yutanoma@dgp.toronto.edu .
