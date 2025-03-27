#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

./build/app/tblite run mol.sdf --method gfn2 --cpcm inf --grad --json 
