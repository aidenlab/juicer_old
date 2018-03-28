#!/bin/bash

res=$1
file=$2
i=$3
topdir=$4
looplist=$5
randomlist=$6

echo $res
echo $file
echo $i
echo $topdir
echo $looplist
echo $randomlist

java -jar hictools.jar dump observed KR $file $i $i BP $res ${topdir}/chr${i}_${res} 

if [ $? -eq 0 ] 
    then  
    awk -v res=$res -f convert_sparse.awk ${topdir}/chr${i}_${res} > ${topdir}/chr${i}_${res}.2
    mv ${topdir}/chr${i}_${res}.2 ${topdir}/chr${i}_${res}
    if [ -n "$looplist" ]
    then
        awk -v res=${res} -v chr=$i -f list_to_matlab.awk "$looplist" > ${topdir}/list${i}_${res}.txt
        awk -v res=${res} -v chr=$i -f list_to_matlab.awk "$randomlist" > ${topdir}/list${i}_${res}_control.txt
        matlab -nodisplay -r "addpath('/Users/muhammadsaadshamim/Desktop/local_arrowhead'); tmp=load('${topdir}/chr${i}_${res}'); tmp=sparse(tmp(:,1), tmp(:,2), tmp(:,3)); list = load('${topdir}/list${i}_${res}.txt'); list1 = load('${topdir}/list${i}_${res}_control.txt'); run_blockbuster(tmp,'${topdir}/chr${i}_${res}_blocks', list, list1); quit;" 
    else
        matlab -nodisplay -r "addpath('/Users/muhammadsaadshamim/Desktop/local_arrowhead'); tmp=load('${topdir}/chr${i}_${res}'); tmp=sparse(tmp(:,1), tmp(:,2), tmp(:,3)); run_blockbuster(tmp,'${topdir}/chr${i}_${res}_blocks'); quit;" 
    fi
    if [ $? -ne 0 ]
    then
        echo "Exit problem in Matlab" 
    else
        if [ -n "$looplist" ]
        then
            rm ${topdir}/list${i}_${res}.txt
            rm ${topdir}/list${i}_${res}_control.txt
            if [ -e ${topdir}/chr${i}_${res}_blocks_scores ]
            then 
                awk -v res=$res -v chr=$i -f convert_list.awk ${topdir}/chr${i}_${res}_blocks_scores > ${topdir}/chr${i}_${res}_scores.txt
                rm ${topdir}/chr${i}_${res}_blocks_scores
            else
                echo "No file ${topdir}/chr${i}_${res}_blocks_scores"
            fi
            if [ -e ${topdir}/chr${i}_${res}_blocks_control_scores ]
            then
                awk -v res=$res -v chr=$i -f convert_list.awk ${topdir}/chr${i}_${res}_blocks_control_scores > ${topdir}/chr${i}_${res}_control_scores.txt
                rm ${topdir}/chr${i}_${res}_blocks_control_scores
            else
                echo "No file ${topdir}/chr${i}_${res}_blocks_control_scores"
            fi
        fi
        if [ -e ${topdir}/chr${i}_${res}_blocks ]
        then
            bash color_blocks.sh ${topdir}/chr${i}_${res}_blocks $i $res > ${topdir}/chr${i}_${res}.txt
            rm ${topdir}/chr${i}_${res}_blocks 
        else
            echo "No file ${topdir}/chr${i}_${res}_blocks"
        fi
        rm ${topdir}/chr${i}_${res}
        exit 0
    fi
else 
    echo "Exit problem in juicebox"
    rm ${topdir}/chr${i}_${res}
    exit 1
fi 
