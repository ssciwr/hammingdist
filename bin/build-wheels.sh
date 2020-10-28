#!/bin/bash
set -e -u -x

# Get a recent version of CMake - the manylinux one is 2.8!
/opt/python/cp38-cp38/bin/pip install cmake
ln -fs /opt/python/cp38-cp38/bin/cmake /usr/bin/cmake

export LD_LIBRARY_PATH=/io/build:$LD_LIBRARY_PATH

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
    fi
}

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" wheel /io/ --no-deps -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
done
