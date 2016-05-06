#!/bin/bash

function get_sources(){
    local l=(
    control.h
    #main_mpi.cc
    general.cc
    general.h
    defs_turb.cc
    defs_turb.h
    funcs.cc
    funcs.h
    stepperbs.cc
    stepperbs.h
    odeintt.h
    odeintt.cc
    )
    echo "${l[@]}"
}

#get_sources
fname_list=$(get_sources)
echo $fname_list
#echo $ll


for fname in ${fname_list[@]}; do
    fname2=$PLA1/src/$fname
    diff_result=$(diff $fname $fname2) # compare text
    #if cmp -s $fname $fname2; then
    if [ "$diff_result" == "" ]; then
        cc="\033[32m" # synced
        ok=1
    else
        cc="\033[31m" # not synced
        ok=0
    fi
    echo -e $cc "--> [$ok] $fname"
done

#echo ${fname_list[0]}
