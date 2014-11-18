#!/bin/bash
# Alignment script. Sets the reference genome and genome ID based on the input
# arguments (default human, DpnII). Optional arguments are the queue for the alignment 
# (default hour), key for menu entry, description, read end, early exit, 
# and the top-level directory (default current directory).
#
# Splits the fastq files, creates jobs to align them, creates merge jobs that 
# wait for the alignment to finish, and creates a final merge and final 
# cleanup job (in case one of the jobs failed).
#
# If all is successful, takes the final merged file, removes name duplicates, 
# removes PCR duplicates, and creates the hic job and stats job.  Final 
# product will be hic file and stats file in the aligned directory.
#
# [topDir]/fastq  - Should contain the fastq files. This code assumes that
#                   there is an "R" in the appropriate files, i.e. *R*.fastq
# From the top-level directory, the following two directories are created:
#
# [topDir]/splits  - Where to write the scratch split files (fastq files and 
#                    intermediate SAM files). This can deleted after execution.
# [topDir]/aligned - Where to write the final output files.
#
# The following globals should be set correctly before proceeding:
#
# splitsize - The number of lines that each split fastq should contain. Larger
#             means fewer files and longer overall, but too small means there
#             are so many jobs that the cluster won't run them
# read1str  - portion of fastq filename that indicates this is the "read 1" 
#             file; used to loop over only the read 1 and within that loop,
#             also align read 2 and merge.  If this is not set correctly, 
#             script will not work. The error will often manifest itself 
#             through a "*" in the name because the wildcard was not able to 
#             match any files with the read1str.
shopt -s extglob

splitsize=6000000 # adjust to match your needs.  40000000 = 1M reads per split
read1str="_R1" # fastq files should look like filename_R1.fastq and filename_R2.fastq
read2str="_R2" # if your fastq files look different, change this value

# unique name for group in case we need to kill jobs
groupname="/a"`date +%s`
# top level directory, can also be set in options
topDir=`pwd`
# default queue, can also be set in options
queue="hour"
# restriction enzyme, can also be set in options
site="DpnII"
# genome ID, default to human, can also be set in options
genomeID="hg19"
# normally both read ends are aligned with long read aligner; if one end is short, this is set
shortreadend=0
# description, default empty
about=""

## Read arguments
usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-q queue] [-s site] [-k key] [-a about] [-R end] [-r] [-e] [-h]"
genomeHelp="   genomeID must be one of \"mm9\" (mouse), \"hg18\" \"hg19\" (human), \"sCerS288c\" (yeast), or \"dMel\" (fly)"
#\n   alternatively, it can be the fasta file of the genome, but the BWA indices must already be created in the same directory"
dirHelp="   [topDir] is the top level directory (default \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
queueHelp="   [queue] is the LSF queue for running alignments (default \"$queue\")"
siteHelp="   [site] must be one of \"HindIII\", \"MseI\", \"NcoI\", \"DpnII\", \"MboI\",\"MspI\", \"HinP1I\", \"StyD4I\", \"SaII\", \"NheI\", \"StyI\", \"XhoI\", \"NlaIII\", or \"merge\" (default \"$site\")"
shortHelp2="   [end]: use the short read aligner on end, must be one of 1 or 2 "
shortHelp="   -r: use the short read version of the aligner (default long read)"
keyHelp="   -k: key for menu item to put this file under"
exitHelp="   -e: early exit; align, sort, merge, and dedup (ends with merged_nodups.txt)"
aboutHelp="   -a: enter description of experiment, enclosed in single quotes"
helpHelp="   -h: print this help and exit"
moreHelp="For more information, type:\n\tgroff -man -Tascii /broad/aidenlab/neva/neva_scripts/align.1 | less"

printHelpAndExit() {
    echo "$usageHelp"
    echo "$genomeHelp"
    echo -e "$dirHelp"
    echo "$queueHelp"
    echo "$siteHelp"
		echo "$keyHelp"
    echo "$shortHelp2"
		echo "$shortHelp"
		echo "$exitHelp"
    echo "$aboutHelp"
    echo "$helpHelp"
    echo -e "$moreHelp"
    exit $1
}

while getopts "d:g:R:k:a:hrefq:s:" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	k) key=$OPTARG ;;
	q) queue=$OPTARG ;;
	s) site=$OPTARG ;;
	R) shortreadend=$OPTARG ;;
	r) shortread=1 ;;  #use short read aligner
	e) earlyexit=1 ;;
  f) fastq=1 ;;
  a) about=$OPTARG ;;
	[?]) printHelpAndExit 1;;
    esac
done


## Set reference sequence based on genome ID
case $genomeID in
    dMel) refSeq="/broad/aidenlab/references/Dmel_release5_first6.fasta";;
    canFam3) refSeq="/seq/references/Canis_lupus_familiaris_assembly3/v0/Canis_lupus_familiaris_assembly3.fasta";;
    mm9) refSeq="/broad/aidenlab/references/Mus_musculus_assembly9_norandom.fasta";;
    hg19) refSeq="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";;
    hg18) refSeq="/seq/references/Homo_sapiens_assembly18/v-1/hg18.fasta";;
		sCerS288c) refSeq="/broad/aidenlab/references/sacCer3.fa";;
    *)  echo "$usageHelp"
        echo "$genomeHelp"
        exit 1
#    *)  refSeq=$genomeID;;
esac

## Check that refSeq exists 
## Note that it really should, given that it is defined above
if [ ! -e "$refSeq" ]; then
    echo "Reference sequence $refSeq does not exist";
    exit 1;
fi


## Set ligation junction based on restriction enzyme
case $site in
    HindIII) ligation="AAGCTAGCTT";;
    MseI)  ligation="TTATAA";;
    DpnII) ligation="GATCGATC";;
    MboI) ligation="GATCGATC";;
    NcoI) ligation="CCATGCATGG";;
    *)  echo "$usageHelp"
				echo "$siteHelp"
				exit 1
esac

## If short read end is set, make sure it is 1 or 2
case $shortreadend in
		0) ;;
		1) ;;
		2) ;;
		*) echo "$usageHelp"
			 echo "$shortHelp2"
			 exit 1
esac

## Hard-coded directory here; in future versions, this should be packaged with script
site_file="/broad/aidenlab/restriction_sites/${genomeID}_${site}.txt"
splitdir=$topDir"/splits"
fastqdir=$topDir"/fastq/*_R*.fastq"
outputdir=$topDir"/aligned"
read1=$splitdir"/*${read1str}*.fastq"
## ARRAY holds the names of the jobs as they are submitted
countjobs=0
declare -a ARRAY

## Check that fastq directory exists and has proper fastq files
if [ ! -d "$topDir/fastq" ]; then
		echo "Directory \"$topDir/fastq\" does not exist."
		echo "Create \"$topDir/$fastq\" and put fastq files to be aligned there."
		echo "Type \"align.sh -h \" for help"
		exit 1
else 
		if stat -t ${fastqdir} >/dev/null 2>&1
				then
				echo "Looking for fastq files...fastq files exist"
		else
				if [ ! -d "$splitdir" ]; then 
						echo "Failed to find any files matching ${fastqdir}"
						echo "Type \"align.sh -h \" for help"
						exit 1			
				fi
		fi
fi

if [ ! -e "$site_file" ]; then
    echo "$site_file does not exist. It must be created before running this script."
    exit 1
fi

## If key option is given, check that the key exists in the menu properties file
if [ ! -z $key ] 
    then
    if ! grep -q -m 1 $key /broad/aidenlab/neva/neva_scripts/hicInternalMenu.properties 
        then 
        echo "Cannot find key $key in hicInternalMenu.properties"
        echo "Please use an existing key or omit the key option"
        exit 1
    fi
fi

## Create output directory
if [ -d "$outputdir" ]; then
    echo "Move or remove directory \"$outputdir\" before proceeding."
		echo "Type \"align.sh -h \" for help"
		exit 1			
fi

mkdir $outputdir

## Create split directory
if [ -d "$splitdir" ]; then
    splitdirexists=1
else
    mkdir $splitdir
fi

## Create temporary directory, used for sort later
if [ ! -d "/broad/hptmp/neva" ]; then
    mkdir /broad/hptmp/neva
    chmod 777 /broad/hptmp/neva
fi

## use LSF cluster on Broad
source /broad/software/scripts/useuse
reuse LSF

echo -e "Aligning files matching $fastqdir\n in queue $queue to genome $genomeID"

## Split fastq files into smaller portions for parallelizing alignment 
## Do this by creating a text script file for the job in "tmp" and then sending it to LSF
if [ ! $splitdirexists ]; then
    echo "Created $splitdir and $outputdir.  Splitting files"
    for i in ${fastqdir}
      do
			echo -e '#!/bin/bash -l' > tmp
			echo -e "#BSUB -q priority" >> tmp
			echo -e "#BSUB -o $topDir/lsf.out\n" >> tmp
			echo -e "#BSUB -g $groupname" >> tmp
			echo -e "#BSUB -J ${groupname}split$i" >> tmp
      filename=$(basename $i)
      filename=${filename%.*}
			echo -e "split -a 3 -l $splitsize -d $i $splitdir/$filename " >> tmp
			ARRAY[countjobs]="${groupname}split${i}"
			countjobs=$(( $countjobs + 1 ))
			bsub < tmp
			rm tmp
    done 


    # Once split succeeds, rename the splits as fastq files
		echo -e '#!/bin/bash -l' > tmp
		echo -e "#BSUB -q priority" >> tmp
		echo -e "#BSUB -o $topDir/lsf.out\n" >> tmp
		echo -e "#BSUB -g $groupname" >> tmp
		echo -e "#BSUB -w \" done(${groupname}split*) \"" >> tmp
		echo -e "#BSUB -J ${groupname}move" >> tmp
    echo -e "for i in $splitdir/*" >> tmp
    echo -e "  do" >> tmp
		echo -e '  mv $i $i.fastq' >> tmp
		echo -e "  done" >> tmp
		bsub < tmp
		ARRAY[countjobs]="${groupname}move"
		countjobs=$(( $countjobs + 1 ))

    # for all the jobs that have been launched so far, create a condition in case
    # any of them exit prematurely
		for (( i=0; i < countjobs; i++ ))
			do
			if [ $i -eq 0 ]; then
					exitjobs="exit(${ARRAY[i]}) "
			else
					exitjobs="$exitjobs || exit(${ARRAY[i]})"
			fi
		done

    # clean up jobs if any fail
		echo -e '#!/bin/bash -l' > tmp
		echo -e "#BSUB -q priority" >> tmp
		echo -e "#BSUB -o $topDir/lsf.out" >> tmp
		echo -e "#BSUB -w \" $exitjobs \"" >> tmp
		echo -e "#BSUB -g ${groupname}kill_splits" >> tmp
		echo -e "#BSUB -J splits_clean_$groupname\n" >> tmp

		echo -e "bkill -g $groupname 0" >> tmp
		bsub < tmp		
		rm tmp
else
    # No need to re-split fastqs if they already exist
    echo -e "Using already created files in $splitdir \n"
fi

## Launch job.  Once split is done, kill the cleanup job, then
## set the parameters for the launch.  Concatenate the script
## "launch.sh" to this script.  This will start the rest of the 
## parallel jobs.
echo -e '#!/bin/bash -l' > tmp2
echo -e "#BSUB -q priority" >> tmp2
echo -e "#BSUB -o $topDir/lsf.out\n" >> tmp2
echo -e "#BSUB -g $groupname" >> tmp2
echo -e "#BSUB -J $groupname$name1$ext" >> tmp2
if [ ! $splitdirexists ]; then
		echo -e "#BSUB -w \"done(${groupname}move)\"" >> tmp2
		echo -e "bkill -g ${groupname}kill_splits 0" >> tmp2
fi
echo "shopt -s extglob" >> tmp2
echo "read1str=\"$read1str\"" >> tmp2
echo "read2str=\"$read2str\"" >> tmp2
echo "read1=\"$read1\"" >> tmp2
echo "queue=\"$queue\"" >> tmp2
echo "topDir=\"$topDir\"" >> tmp2
echo "groupname=\"$groupname\"" >> tmp2
echo "refSeq=\"$refSeq\"" >> tmp2
echo "site_file=\"$site_file\"" >> tmp2
echo "splitdir=\"$splitdir\"" >> tmp2
echo "outputdir=\"$outputdir\"" >> tmp2
echo "ligation=\"$ligation\"" >> tmp2
echo "genomeID=\"$genomeID\"" >> tmp2
echo "about=\"$about\"" >> tmp2
echo "shortreadend=$shortreadend" >> tmp2
flags="-g $genomeID -d $topDir -s $site -q $queue -a '$about'"
echo "flags=\"$flags\" " >> tmp2
if [ $shortread ]; then
		echo "shortread=$shortread" >> tmp2
		echo "flags=\"$flags -r \"" >> tmp2
fi
if [ $shortreadend -gt 0 ]; then
		echo "flags=\"$flags -R $shortreadend \"" >> tmp2
fi
if [ $earlyexit ]; then
		echo "earlyexit=$earlyexit" >> tmp2
		echo "flags=\"$flags -e \"" >> tmp2
fi
if [ ! -z $key ]; then
		echo "key=\"$key\"" >> tmp2
		echo "flags=\"$flags -k $key \"" >> tmp2
fi

if [ $splitdirexists ]; then
		echo "splitdirexists=$splitdirexists" >> tmp2
fi

cat tmp2 /broad/aidenlab/neva/neva_scripts/launch.sh > tmp3
echo "Starting job to launch other jobs once splitting is complete"
bsub < tmp3
rm tmp2 tmp3
