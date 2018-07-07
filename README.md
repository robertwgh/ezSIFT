## ezSIFT: An easy-to-use stanalone SIFT library 

[![License][license-img]][license-url] [![Build Status](https://travis-ci.com/robertwgh/ezSIFT.svg?branch=master)](https://travis-ci.com/robertwgh/ezSIFT)

The SIFT (scale-invariant feature transform) algorithm is considered to be one of the most robust local feature detector and description methods. Most of the open-source SIFT implementations rely on some 3rd-party libraries. Some of them even rely on a few different large libraries. These dependencies make the installation, compilation and usage not easy.

The ezSIFT library provides a standalone and lightweight SIFT implementation written in C/C++. The ezSIFT is self-contained, and does not require any other libraries. So it is easy to use and modify. Besides, the implementation of the ezSIFT is straightforward and easy to read. 

The ezSIFT project was original hosted at SourceForge ([link](https://sourceforge.net/projects/ezsift)). We recently migrated it to github to better maintain the code development. 

### Documentation
Please read [ezSIFT Wiki page](https://github.com/robertwgh/ezSIFT/wiki) for details.

### Examples
I also provide two examples showing how to use this library:

* `examples/feature_extract`: detect keypoints and extract feature descriptor from a single image.
* `examples/image_match`: detect keypoints and extract features from two images and perform feature matching. 

### How to build
#### Mac OS
Follow the following instructions:
```Bash
cd platforms/desktop/
mkdir build
cd build
cmake ..
make
```
Then you can find the built binary under `build/bin` directory. Run the two demos like this:

```bash
./image_match img1.pgm img2.pgm
./feature_extract img1.pgm
```

Or, you can use the following instruction to generate Xcode project:
```Bash
cd platforms/desktop/
mkdir build
cd build
cmake .. -GXcode
```
Then, open `ezsift.xcodeproj` project to build.

#### Windows (Visual Studio)
```Bash
cd platforms/desktop/
mkdir build
cd build
cmake ..
```
Upon finishing the above instructions, A Visual Studio solution `ezsift.sln` will be generated. Open the solution file and build all projects, the static lib and examples will be built. 

#### Android NDK native
1. Please install the latest Android NDK package, and add NDK root folder to your system environment PATH to expose `ndk-build` tool to the shell. 
2. Go to `platforms/android`.
3. run `./build.sh`
4. After build is finished, you will find the binaries under `build/android/bin/`.
5. Connect an Android device to your computer using ADB. Use `install_and_run.sh` to install the binaries and run them on the device.

### License

    Copyright 2013 Guohui Wang

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.


[license-url]: https://github.com/robertwgh/ezSIFT/blob/master/LICENSE
[license-img]: https://img.shields.io/badge/License-Apache%202.0-blue.svg
