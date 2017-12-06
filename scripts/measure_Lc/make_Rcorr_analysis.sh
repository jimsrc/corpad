#!/bin/bash
#++++++++++++++++++++++++++++++++++++++++++++++
ok_run=0        # default: don't run
ok_fit=0        # default: don't fit lambda_c value
#--- parse args
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -odir|--odir)
        dir_dst="$2"
        echo -e "\n [*] Output dir: ${dir_dst}\n"
        shift # past argument
        shift # past value
        ;;
        -run|--run)
        ok_run=1
        echo -e "\n [*] we'll run dB model to build R(r,th) function!\n"
        shift # past arg
        ;;
        -fit|--fit) # fit lambda_c value (for each run)
        ok_fit=1
        echo -e "\n [*] we'll fit correlation lengths (lambda_c) values\n"
        shift # past arg
        ;;
        *)    # unknown option
        echo -e "\n [-] unknown argument: $key\n"
        exit 1
        ;;
    esac
done
#++++++++++++++++++++++++++++++++++++++++++++++

# number of dB realizations
Nrlz=256

# NOTE: reasonables values for Lc are from <=0.1. This is because
# a value of, for example, 0.3 would be too close to 1AU=lmax, which
# is the largest wavelength of the modeled turbulence.
# In fact, for cases with Lc>=0.3 we obtain weird correlation function.
Lc_values=(0.1 0.03 0.01 0.00562 0.003)

#--- build R(r,th)
for Lc in ${Lc_values[@]}; do
    fout=${dir_dst}/Rbidim_Lc.$(printf "%1.2e" $Lc).h5
    echo " +++++++++++++++++++++++++++++++++++++++"
    [[ ${ok_run} -eq 1 ]] && {
        # ok, we'll run
        echo -e " [*] output: $fout\n":
        ./test.py --pdb -- \
            --bidim \
            -o $fout \
            -sigma 0.3 \
            -slab 0.2 \
            -Lc $Lc \
            -Nrlz $Nrlz \
            -nLc 80 \
            -nth 31 \
            -ndr 64 \
            -ndrLc 3.5 \
            2> /dev/null \
            || {
                echo -e " [-] sthing wen't wrong!";
                exit 1;
            };
    }

    # generate isotrpized R(r,th)-->R(r), then fit the value of the
    # correlation length lambda_c && plot
    ffig=${dir_dst}/Rbidim_Lc.$(printf "%1.2e" $Lc)_Nrlz.$Nrlz.png
    [[ ${ok_fit} -eq 1 ]] && {
        echo -e "\n [*] generating isotropized version of R(r,th) && fit lambda_c ...\n";
        ./plot_R.py --pdb -- \
            -i $fout \
            -nlevels 16 \
            -fig $ffig \
            || exit 1;
    }
done

#EOF
