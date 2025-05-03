# GeometryTemplate
A template project to get started with geometry-central and Polyscope. Copied from Nick's template [here](https://github.com/nmwsharp/gc-polyscope-project-template), but with code for googletest added in.

This repo is set up as a template, but it seems the submodules don't get copied to the new project. Run `setup.sh` to set up the required submodules and to run cmake. Setup also sets up cmake to export compile commands for `clang-format`. To build the code, you can run
```
cd build
make -j7
```

Then run the code with
```
bin/run /path/to/a/mesh
```

Run the tests with
```
bin/test
```
