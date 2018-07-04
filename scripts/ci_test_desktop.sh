# This script will be called by CI from the root directory
echo "run examples for desktop"
cd ./platforms/desktop/build/bin
./image_match img1.pgm img2.pgm
./feature_extract img1.pgm