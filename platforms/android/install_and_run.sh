#!/bin/sh
adb shell 'rm -rf /data/local/tmp/ezsift'
adb shell "mkdir -p /data/local/tmp/ezsift/"
adb push bin /data/local/tmp/ezsift/
adb push ../../data/img1.pgm /data/local/tmp/ezsift/bin/
adb push ../../data/img2.pgm /data/local/tmp/ezsift/bin/
adb shell 'chmod 777 -R /data/local/tmp/ezsift/bin'

adb shell "adb shell "cd /data/local/tmp/ezsift/bin && \
           ./feature_extract img1.pgm""

adb shell "adb shell "cd /data/local/tmp/ezsift/bin && \
	       ./image_match img1.pgm img2.pgm""