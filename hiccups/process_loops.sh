#!/bin/bash

# set genome to human by default
genomeID="hg19"
# set looplist to 0 by default
pixellist="0"
# set radius to 20000 by default
radius=20000
# set kr input file to 0 by default
kr_input="0"
# set res to 10000 by default
res=10000
# set fdrsum to 0.02 by default
fdrsum="0.02"
# set output to output by default
outputdir="output"
# set color to cyan by default
color="0,255,255"
# set color switch to no by default
colorswitch="n"
# set o/e threshold1 to 1.5 by default
oethreshold1="1.5"
# set o/e threshold2 to 1.75 by default
oethreshold2="1.75"
# set o/e threshold3 to 2 by default
oethreshold3="2"


## Read arguments
usageHelp="Usage: ${0##*/} -g genomeID [-p pixellist] [-r radius] [-s resolution] [-f fdrsum] [-k kr_input] [-o output] [-c color] [-w colorswitch] [-n name] [-h]"
genomeHelp="   genomeID must be one of \"mm9\" (mouse) or \"hg19\" (human); default  is \"$genomeID\""
#\n   alternatively, it can be the fasta file of the genome, but the BWA indices must already be created in the same directory"
pixelHelp="   [pixellist] must be a list of pixels returned by HiCCUPS"
radiusHelp="   [radius] must be an integer (initial radius size for clustering); default is $radius"
resolutionHelp="   [resolution] should be either 5000, 10000, or 25000; default is $res"
fdrsumHelp="   [fdrsum] should be a float, the max value that all four fdr values can sum to for a singleton cluster; default is $fdrsum"
krinputHelp="   [kr_input] should be a file containing the chromosome name in field 1 and the path to the kr norm vectors in field2"
outputHelp="   [output] should be a directory to put all output files in"
colorHelp="   [color] should be an rgb string for peaks in appear in the viewer; default is $color"
colorswitchHelp="   [colorswitch] should be y or n; if y, color for peaks in the viewer will vary based on fdr"
nameHelp="   [name] should be the name for the final peak file"
helpHelp="   -h: print this help and exit"

printHelpAndExit() {
    echo "$usageHelp"
    echo "$genomeHelp"
    echo "$pixelHelp"
    echo "$radiusHelp"
    echo "$resolutionHelp"
    echo "$fdrsumHelp"
    echo "$krinputHelp"
    echo "$outputHelp"
    echo "$colorHelp"
    echo "$colorswitchHelp"
    echo "$nameHelp"
    echo "$helpHelp"
    exit $1
}

while getopts "g:p:hr:s:f:k:o:c:w:n:" opt; do
    case $opt in
        g) genomeID=$OPTARG ;;
        h) printHelpAndExit 0;;
        p) pixellist=$OPTARG ;;
        r) radius=$OPTARG ;;
        s) res=$OPTARG ;;
        k) kr_input=$OPTARG ;;
        f) fdrsum=$OPTARG ;;
	o) outputdir=$OPTARG ;;
	c) color=$OPTARG ;;
	w) colorswitch=$OPTARG ;;
	n) finalname=$OPTARG ;;
        [?]) printHelpAndExit 1;;
    esac
done

if [ $pixellist == 0 ]
then
        echo "You must input a pixellist"
        exit $1
fi

if [ $kr_input == 0 ]
then
        echo "You must input a kr input file"
        exit $1
fi

mkdir -p "$outputdir"
/raid10/suhas/remove_lowmapq.py "$kr_input" "$pixellist" "$res" > "$outputdir""/loops_10fdr_mapq_corrected.txt"
if [ $genomeID == "hg19" ]
then
	awk '{if ($1=="X") {print 23 "\t" $2 "\t" 23 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17} else {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17}}' "$outputdir""/loops_10fdr_mapq_corrected.txt" > "$outputdir""/loops_10fdr_mapq_corrected_23X.txt"
	/raid10/suhas/unique_peaks_centroid.py "$outputdir""/loops_10fdr_mapq_corrected_23X.txt" "$radius" "$res" "$color" "$colorswitch" > "$outputdir""/loops_10fdr_mapq_corrected_unique.txt"
	awk -v t1=$oethreshold1 -v t2=$oethreshold2 -v t3=$oethreshold3 -v f=$fdrsum 'BEGIN{print "chr1" "\t" "x1" "\t" "x2" "\t" "chr2" "\t" "y1" "\t" "y2" "\t" "color" "\t" "o" "\t" "e_bl" "\t" "e_donut" "\t" "e_h" "\t" "e_v" "\t" "fdr_bl" "\t" "fdr_donut" "\t" "fdr_h" "\t" "fdr_v" "\t" "num_collapsed" "\t" "centroid1" "\t" "centroid2" "\t" "radius"}{if ($8>t2*$9&&$8>t2*$10&&$8>t1*$11&&$8>$12*t1&&($8>t3*$9||$8>t3*$10)&&($21>1||($17+$18+$19+$20)<=f)) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24}}' "$outputdir""/loops_10fdr_mapq_corrected_unique.txt" > "$outputdir""/""$finalname"
elif [ $genomeID == "mm9" ]
then
	awk '{if ($1=="X") {print 23 "\t" $2 "\t" 23 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17} else {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17}}' "$outputdir""/loops_10fdr_mapq_corrected.txt" > "$outputdir""/loops_10fdr_mapq_corrected_23X.txt"
	/raid10/suhas/unique_peaks_centroid_mouse.py "$outputdir""/loops_10fdr_mapq_corrected_23X.txt" "$radius" "$res" "$color" "$colorswitch" > "$outputdir""/loops_10fdr_mapq_corrected_unique.txt"
	awk -v t1=$oethreshold1 -v t2=$oethreshold2 -v t3=$oethreshold3 -v f=$fdrsum 'BEGIN{print "chr1" "\t" "x1" "\t" "x2" "\t" "chr2" "\t" "y1" "\t" "y2" "\t" "color" "\t" "o" "\t" "e_bl" "\t" "e_donut" "\t" "e_h" "\t" "e_v" "\t" "fdr_bl" "\t" "fdr_donut" "\t" "fdr_h" "\t" "fdr_v" "\t" "num_collapsed" "\t" "centroid1" "\t" "centroid2" "\t" "radius"}{if ($8>t2*$9&&$8>t2*$10&&$8>t1*$11&&$8>$12*t1&&($8>t3*$9||$8>t3*$10)&&($21>1||($17+$18+$19+$20)<=f)) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24}}' "$outputdir""/loops_10fdr_mapq_corrected_unique.txt" > "$outputdir""/""$finalname"
elif [ $genomeID == "galGal4" ]
then
	awk '{if ($1=="Z") {print 29 "\t" $2 "\t" 29 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17} else {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17}}' "$outputdir""/loops_10fdr_mapq_corrected.txt" > "$outputdir""/loops_10fdr_mapq_corrected_29Z.txt"
	/raid10/suhas/unique_peaks_centroid_chicken.py "$outputdir""/loops_10fdr_mapq_corrected_29Z.txt" "$radius" "$res" "$color" "$colorswitch" > "$outputdir""/loops_10fdr_mapq_corrected_unique.txt"
	awk -v t1=$oethreshold1 -v t2=$oethreshold2 -v t3=$oethreshold3 -v f=$fdrsum 'BEGIN{print "chr1" "\t" "x1" "\t" "x2" "\t" "chr2" "\t" "y1" "\t" "y2" "\t" "color" "\t" "o" "\t" "e_bl" "\t" "e_donut" "\t" "e_h" "\t" "e_v" "\t" "fdr_bl" "\t" "fdr_donut" "\t" "fdr_h" "\t" "fdr_v" "\t" "num_collapsed" "\t" "centroid1" "\t" "centroid2" "\t" "radius"}{if ($8>t2*$9&&$8>t2*$10&&$8>t1*$11&&$8>$12*t1&&($8>t3*$9||$8>t3*$10)&&($21>1||($17+$18+$19+$20)<=f)) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24}}' "$outputdir""/loops_10fdr_mapq_corrected_unique.txt" > "$outputdir""/""$finalname"
fi
