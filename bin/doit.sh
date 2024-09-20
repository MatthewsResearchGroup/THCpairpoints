#!/bin/bash


# modify the method and basis set as you want
METHOD="MP3 CCSD MP2" 
BASIS="PVDZ AUG-PVDZ PVTZ REDUCED_PVTZ"

ZMATfile=./ZMAT
dummyfile=./rundummy

# add the loops for molecules in the further
for calctype in $METHOD
do
    for basisset in $BASIS
    do  
        jobname=$basisset-$calctype
        mkdir -p $jobname
        cd $jobname
        cp ../$ZMATfile  ./
        cp ../$dummyfile ./runjob.sh
        # replace some keywords in dummyfile
        sed -i s/dummyrun/$jobname/ runjob.sh
        sed -i s/dummyworkdir/'$PWD'/ runjob.sh
        if [ "$calctype" = "MP2" ]; then
            sed -i s/dummyextract// runjob.sh
            sed -i s/aaa/LS/ ZMAT
        else
            sed -i s/dummyextract/'python3 $EXTRACT'/ runjob.sh
            sed -i s/aaa/OFF/ ZMAT
        fi
        
        # replace some keywords in ZMAT
        sed -i s/xxx/$calctype/ ZMAT
        # since the reduced-cc-PVTZ is a little bit different, we need to replace with the string "PVTZ,GENBAS_1=320,GENBAS_2=4320"
        if [ "$basisset" = "REDUCED_PVTZ" ]; then
            sed -i s/yyy/PVTZ,GENBAS_1=320,GENBAS_2=4320/ ZMAT
        else
            sed -i s/yyy/$basisset/ ZMAT
        fi 

        cd ../ 
    done
done

