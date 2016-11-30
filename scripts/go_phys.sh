#!/bin/bash
F1=../out/h.001/post.h5
F2=../out/h.002/post.h5
F3=../out/h.004/post.h5

./make_phys.py --pdb -- --LcSlab 0.01  --xi 1.0  --Bo 5e-5  --fnames $F1 $F2 $3

#EOF
