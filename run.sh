#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

./build/app/tblite run coord.xyz --no-restart --method gfn2 --lpb 7.0 --kappa 10 
