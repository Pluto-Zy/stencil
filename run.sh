#!/bin/bash

bsub -I -b -q q_ustc -N 1 -np 1 -cgsp 64 -cache_size 0 -PARSE slave ./build/bin/stencil_main -m 400 -b 50 -i 1000 -w 1