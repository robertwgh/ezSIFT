##ezSIFT: An easy-to-use stanalone SIFT library.
URL: https://sourceforge.net/projects/ezsift/

###Version
Created on 9/16/2013.
Modified on 3/27/2015.

###Examples
I also provide two examples showing how to use this library:
`examples/feature_extract.cpp`: detect keypoints and extract feature descriptor from a single image.
`examples/image_match.cpp`: detect keypoints and extract features from two images and perform feature matching. 

###How to build
####For Visual Studio 2010/2012
1. Go to `build/vs2010` or `build/vs2012`, open `ezsift.sl`n. 
2. Build solution. 
3. You will find executable files under `build/vs2010/bin` or `build/vs2012/bin`.

####For Android NDK native mode
1. Please install Android NDK package, and add NDK root folder to your system environment PATH. 
2. Go to `build/android`.
3. run `./build.sh`
4. You will find the binaries under `build/android/bin/`.
