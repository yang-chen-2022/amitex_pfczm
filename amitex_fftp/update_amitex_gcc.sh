#!/bin/bash

git checkout install

git pull

sed -i "/^#FC\=/c\FC=gfortran" install

./install

cd validation

qsub script_tests_marquises

cd ..

more validation/tests.log
more validation/time.log
