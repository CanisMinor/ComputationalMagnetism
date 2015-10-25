#!/bin/bash
rm magnn.dat
rm CMCFePt
rm magn-vs-T.dat
rm neighbours.dat
g++ RotMat.cpp NeighbReadIn.cpp ExchangeReadIn.cpp cmoncar.cpp -o CMCFePt  -O3 -Wall -pedantic -msse3 

#ifort cmoncar.f90 -o heun -O3 -static -ipo -axT 
./CMCFePt
#time ./CMCFePt	
