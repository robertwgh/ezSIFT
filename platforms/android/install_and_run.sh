#!/bin/sh
adb shell 'rm -rf /data/local/tmp/ezsift'
adb shell "mkdir -p /data/local/tmp/ezsift/"
adb push bin /data/local/tmp/ezsift/
adb push ../../data/img1.pgm /data/local/tmp/ezsift/bin/
adb push ../../data/img2.pgm /data/local/tmp/ezsift/bin/
adb shell 'chmod 777 -R /data/local/tmp/ezsift/bin'

echo "\n******************************"
echo "* Running feature extraction *"
echo "******************************"
adb shell "cd /data/local/tmp/ezsift/bin && \
           export LD_LIBRARY_PATH=./:$LD_LIBRARY_PATH && \
           ./feature_extract img1.pgm"

echo "\n***********************"
echo "* Running image match *"
echo "***********************"
adb shell "cd /data/local/tmp/ezsift/bin && \
           export LD_LIBRARY_PATH=./:$LD_LIBRARY_PATH && \
           ./image_match img1.pgm img2.pgm"