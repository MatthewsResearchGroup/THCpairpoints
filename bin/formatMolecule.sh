#!/bin/bash
# post processing for molecules grid

SECONDS=0
mkdir -p data

types=(adz . tz rtz)
dirs=(AUG-PVDZ-CCSD AUG-PVDZ-MP2 AUG-PVDZ-MP3 PVDZ-CCSD PVDZ-MP2 PVDZ-MP3 PVTZ-CCSD PVTZ-MP2 PVTZ-MP3 REDUCED_PVTZ-CCSD REDUCED_PVTZ-MP2 REDUCED_PVTZ-MP3)

echo "Starting copy files..."
for i in {0..3}
do
    cd data
    type=${types[$i]}
    exit="./.."
    if [ "$type" != "." ]; then
        mkdir -p $type
        cd $type
        exit="./../.."
    fi
    mkdir -p grid
    mkdir -p ccsd
    
    # copy MP3 data
    cp $exit/${dirs[$i*3+2]}/fa.dat ./
    cp $exit/${dirs[$i*3+2]}/fi.dat ./
    cp $exit/${dirs[$i*3+2]}/v.dat ./
    cp $exit/${dirs[$i*3+2]}/t.dat ./
    
    # change to grid directory and copy some files
    cd grid
    cp $exit/../${dirs[$i*3+1]}/XA.dat ./
    cp $exit/../${dirs[$i*3+1]}/XI.dat ./
    cp $exit/../${dirs[$i*3+1]}/grid.dat ./
    cp $exit/../${dirs[$i*3+1]}/pvt.dat ./

    # change to ccsd directory and copy realted files
    cd ../ccsd
    cp $exit/../${dirs[$i*3+0]}/fa.dat ./
    cp $exit/../${dirs[$i*3+0]}/fi.dat ./
    cp $exit/../${dirs[$i*3+0]}/v.dat ./
    cp $exit/../${dirs[$i*3+0]}/t.dat ./
    
    cd $exit/..
    echo $PWD
done
echo "Copy done!"
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
