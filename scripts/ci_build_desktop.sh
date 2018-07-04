echo "build for desktop"
cd $(pwd)/../platforms/desktop
mkdir -p build
cd build
cmake ..
make -j4