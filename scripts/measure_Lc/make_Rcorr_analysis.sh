#!/bin/bash

dir_dst="./tmp"  # "."
Nrlz=256
#Lc_values=(3.0 1.0 0.3 0.1 0.03 0.01)
Lc=1.0

#--- build R(r,th)
#for Lc in ${Lc_values[@]}; do
fout=${dir_dst}/Rbidim_Lc.$(printf "%1.2f" $Lc).h5
echo " +++++++++++++++++++++++++++++++++++++++"
echo " [*] output: $fout\n"
./test.py -- \
    --bidim \
    -o $fout \
    -sigma 0.3 \
    -slab 0.2 \
    -Lc $Lc \
    -Nrlz $Nrlz \
    -nLc 80 \
    -nth 31 \
    -ndr 32 \
    -ndrLc 3.5 \
    2> /dev/null \
    || {
        echo -e " [-] sthing wen't wrong!";
        exit 1;
    }

# make fig && fit
ffig=${dir_dst}/Rbidim_Lc.$(printf "%1.2f" $Lc)_Nrlz.$Nrlz.png
./plot_R.py --pdb -- \
    -i $fout \
    -nlevels 16 \
    -fig $ffig \
    || exit 1
#done

##--- generate isotropized version of R(r,th)-->R(r)
#echo    " +++++++++++++++++++++++++++++++++++++++++"
#echo -e "\n [*] generating isotropized versions..."
#for Lc in ${Lc_values[@]}; do
#    finp=${dir_dst}/Rbidim_Lc.$(printf "%1.2f" $Lc).h5
#    ffig=${dir_dst}/Rbidim_Lc.$(printf "%1.2f" $Lc).png
#    echo    " +++++++++++++++++++++++++++++++++++++++++"
#    echo " [*] $finp "
#    ./plot_R.py --pdb -- \
#        -i $finp \
#        -fig $ffig \
#        || exit 1
#done

#EOF
