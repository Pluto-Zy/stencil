#!/bin/bash

bsub -I -b -q q_ustc -N 1 -np 1 -cgsp 64 -cache_size 0 ./build/bin/stencil_main -s 400 -b 50 -i 1000 -r 1 -R 3 -m DMA RMAWithoutSync -c