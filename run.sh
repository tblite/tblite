#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

./build/app/tblite run ibu.xyz --lpb water --kappa 5.0 
