#!/usr/bin/env bash

rm -rf build bin lib
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j12