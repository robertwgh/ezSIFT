#!/bin/sh
ndk-build
mkdir -p bin/
cp -rf libs/armeabi-v7a/* bin/