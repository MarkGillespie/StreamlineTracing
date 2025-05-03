#!/bin/bash

mkdir deps
cd deps
git submodule add -b master git@github.com:MarkGillespie/polyscope.git
git submodule add -b master git@github.com:MarkGillespie/geometry-central.git
git submodule add -b v1.12.x git@github.com:google/googletest
git submodule update --init --recursive
cd ..
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
