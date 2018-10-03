#!/bin/bash


#install OGDF
cd OGDF
bash makeMakefile.sh
make -j 16
cd ..

#install METACARVEL
make
