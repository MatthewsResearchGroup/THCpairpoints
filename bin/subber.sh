#!/bin/bash
for dir in */; do 
  cd $dir
  sbatch runjob.sh 
  cd ..
done
