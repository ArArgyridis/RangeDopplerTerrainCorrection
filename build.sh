#! /bin/sh

[[ -d build ]] || mkdir build
[[ -d bin ]] || mkdir bin

cd build

cmake ../

make 

mv RangeDopplerTerrainCorrection ../bin/
