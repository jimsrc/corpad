#!/bin/bash
prefix=h_lmax #h_lmax #h
F1=../out/${prefix}.001/post.h5
F2=../out/${prefix}.002/post.h5
F3=../out/${prefix}.004/post.h5
F4=../out/${prefix}.005/post.h5
#F4=../out/${prefix}.004/post.h5
#F5=../out/${prefix}.005/post.h5

./make_phys.py --pdb -- --xi 1.0  \
    --fnames $F1 $F2 $F3 $F4 \
    --ro .2 .3 .7 .9 \
    --figname ${prefix}_test.png

#EOF
