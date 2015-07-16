#!/bin/bash

mpiCC create_mesh.cpp -fopenmp -o create_mesh -I/work/klicpera/simmodeler/include -L/work/klicpera/simmodeler/lib/x64_rhel5_gcc41 -lSimMeshing -lSimDiscrete -lSimMeshTools -lSimModel
