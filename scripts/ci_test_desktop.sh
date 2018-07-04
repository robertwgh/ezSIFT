echo "run examples for desktop"
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${dir}/../platforms/desktop/build/bin
./image_match img1.pgm img2.pgm
./feature_extract img1.pgm