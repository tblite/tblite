#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

./build/app/tblite run mindless.xyz --iterations 1000 --cosmo 78 --grad --json 
