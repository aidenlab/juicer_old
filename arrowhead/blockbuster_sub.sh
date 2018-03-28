#!/bin/bash

res=$1
file=$2
i=$3;

java -jar hictools.jar dump observed KR $file $i $i BP $res > chr${i}_${res} 
if [ $? -eq 0 ] 
    then  
    awk -v res=$res -f convert_sparse.awk chr${i}_${res} > chr${i}_${res}.2
    mv chr${i}_${res}.2 chr${i}_${res}
    matlab -nodisplay -r "addpath('/Users/muhammadsaadshamim/Desktop/local_arrowhead'); tmp=load('chr${i}_${res}'); tmp=sparse(tmp(:,1), tmp(:,2), tmp(:,3)); run_blockbuster(tmp,'chr${i}_${res}_blocks'); quit;" 
    if [ $? -ne 0 ]
    then
        echo "Exit problem in Matlab or color blocks" 
        exit 1
    else  
        color_blocks.sh chr${i}_${res}_blocks $i $res > chr${i}_${res}.txt
        rm chr${i}_${res}_blocks 
    fi
else 
    echo "Exit problem in hictools"
fi 
