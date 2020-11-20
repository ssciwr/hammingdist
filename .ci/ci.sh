#!/usr/bin/env bash

# compile code & run tests

pwd
ls
echo "$BUILD_TYPE"
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j2
make test
cd ..

