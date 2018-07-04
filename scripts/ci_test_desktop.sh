echo "run examples for desktop"
cd $(pwd)/../platforms/desktop/build/bin
./image_match img1.pgm img2.pgm
./feature_extract img1.pgm