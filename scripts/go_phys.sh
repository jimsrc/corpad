#!/bin/bash
F1=../out/h.001/post.h5
F2=../out/h.002/post.h5
F3=../out/h.003/post.h5
F4=../out/h.004/post.h5
F5=../out/h.005/post.h5

./make_phys.py --pdb -- --xi 1.0  \
    --fnames $F1 $F2 $F3 $F4 $F5 \
    --ro .2 .3 .5 .7 .9

#EOF
