echo "build for desktop"
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${dir}/../platforms/desktop
mkdir -p build
cd build
cmake ..
make -j4