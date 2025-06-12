# Make sure to specify $prefix in your shell environment before running this script
# $prefix specifies where you want Tasmanian installed
# CMAKE_CUDA_COMPILER: The location of your CUDA nvcc
# MAKE_CUDA_ARHITECTURES: This should cover many common architectures. Check your own architecture to make sure
# DBLA_VENDOR: Make sure that you specify your BLAS vendor

cmake -DCMAKE_INSTALL_PREFIX=$prefix \
    -DCMAKE_CUDA_COMPILER=/your/path/to/cuda/version/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES="60;70;75;89;90"\
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DTasmanian_ENABLE_RECOMMENDED=ON \
    -DTasmanian_ENABLE_PYTHON=OFF \
    -DBLA_VENDOR=OpenBLAS\
    -DDEBUG_FIND=ON\
    -DTasmanian_ENABLE_CUDA=ON\
    path/to/tasmanian/source

# Full pipeline (from the original Tasmanian repo)
## export prefix="path/to/desired/install/location"
##  mkdir Build
##  cd Build
##  cmake <options> <path-to-Tasmanian-source>
##  make
##  make test
##  make install
##  make test_install
