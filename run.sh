#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

./build/app/tblite run coord.xyz --no-restart --method gfn2 --cosmo 7.0 
