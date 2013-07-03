#!/bin/bash

FILE="final1"

(cd .. && make)
rm $FILE.exr
./pbrt scenes/$FILE.pbrt
exrdisplay $FILE.exr


