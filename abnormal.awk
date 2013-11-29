# Takes a SAM file, looks for interesting abnormal chimeric reads.
# only those reads mapped to chromosomes 1-24 with MAPQ >= 10, and
# with at least a triple where the minimum distance >= 20Kb
#
function abs(v) {
    return v<0?-v:v;
}
BEGIN{
  OFS="\t";
  mapq_thresh=-1;
  tottot = -1; # will count first non-group
}
{
  # input file is sorted by read name.  Look at read name to group 
  # appropriately
  split($1,a,"/");
  if(a[1]==prev){
    # move on to next record.  look below this block for actions that occur
    # regardless.
    count++;
  }
  else {
    # deal with read pair group
    tottot++;
#    printme = 0;
    printme = 1;
    j = 0;
    for (i in c) {
      # line contains ligation junction
#      if (c[i] ~ /AAGCTAGCTT/) {
#	printme = 1;
#      }
      split(c[i], tmp);
      mapped[j] = tmp[3] ~ /^[1-9,X,Y,M][0-9,T]?$/ || tmp[3] ~ /^chr[1-9,X,Y,M][0-9,T]?$/;
      name[j] = tmp[1];
      str[j] = tmp[2];
      chr[j] = tmp[3];
      pos[j] = tmp[4];
      m[j] = tmp[5];
      cig[j] = tmp[6];
      seq[j] = tmp[10];
      j++;
    }
    len = j;
    bigdist = 0;
    for (j=0; j<len; j++) {
        printme = printme && mapped[j] && (m[j] >= 10);
        for (k=j+1; k<len; k++) {
            if (abs(pos[j]-pos[k]) >= 20000) {
                bigdist++;
            }
        }
    }
    if (printme && bigdist >= 3) {
      for (j = 0; j < len; j++) {
          printf("%s %d %s %d %d %s %s ", name[j], str[j], chr[j], pos[j], m[j], cig[j], seq[j]);
      }
      printf("\n");
    }

    # reset variables
    delete c;
    count=1;
    prev=a[1];
  }
  # these happen no matter what, after the above processing
  c[count] = $0;
}
