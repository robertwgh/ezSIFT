#!/bin/sh
ndk-build
mkdir -p bin/armeabi-v7a/
cp -rf libs/armeabi-v7a/* bin/armeabi-v7a/