BEGIN{
    FS=",";
    OFS="\t";
    print "chr1","x1","x2","chr2","y1","y2","block score";
}
{
    $1=$1-1;
    $2=$2-1;
    $3=$3-1;
    $4=$4-1;
    $1=$1*res;
    $2=$2*res;
    $3=$3*res;
    $4=$4*res;
    print chr,$1,$2,chr,$3,$4,$5;
}
