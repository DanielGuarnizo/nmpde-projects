#!/bin/bash

# Delete build folder if it exists
if [ -d build ]; then
    rm -rf build
fi
# Load neccessary library
module load gcc-glibc dealii

# Create build directory and compile
mkdir build && cd build
cmake ..
cmake --build .

# Run executables
# ./neuro_disease_3D
./neuro_disease_1D