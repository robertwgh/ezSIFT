# This script will be called by CI from the root directory
echo "build for desktop"
cd ./platforms/desktop
mkdir -p build
cd build
cmake ..
make -j4