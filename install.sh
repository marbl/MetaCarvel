#!/bin/bash

set -euo pipefail

#install OGDF
cd OGDF
./makeMakefile.py
make -j 16
cd ..

#install METACARVEL
make
