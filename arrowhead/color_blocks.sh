#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Usage: color_blocks.sh file chr res";
    exit 1;
fi

awk -v res=$3 -v chr=$2 'BEGIN{
    OFS="\t"; 
    FS=",";
    red=255;
    green=255;
    blue=0;
}
NR==1 {
    printf "chr1\tx1\tx2\tchr2\ty1\ty2\tcolor";
    for (i=3; i<=NF; i++) {
       printf "\tf%d",  i-2;
    }
    printf "\n";
}
{
    color=sprintf("%d,%d,%d",red,green,blue);  
    printf chr"\t"($1-1)*res"\t"$2*res"\t"chr"\t"($1-1)*res"\t"$2*res"\t"color;
    for (i=3; i<=NF; i++) {
       printf "\t%s",  $i;
    }
    printf "\n";
}' $1
if [ $? -ne 0 ]; then
    exit 1;
fi
