#!/bin/bash

set -e

# Remove all sorts of old build artifacts
rm -rf build dist wheelhouse hammingdist.egg-info

export DOCKER_IMAGE=quay.io/pypa/manylinux2010_x86_64
export PLAT=manylinux2010_x86_64

docker pull $DOCKER_IMAGE
docker run --rm -e PLAT=$PLAT -v `pwd`:/io $DOCKER_IMAGE /io/bin/build-wheels.sh

python -m pip install --upgrade twine
twine upload --repository testpypi wheelhouse/*
