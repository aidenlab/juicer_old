#!/usr/bin/env python

import sys
import numpy as np

kr_factors={}

input_file=open(sys.argv[1],'r')
input_line=input_file.readline()

while input_line!="":
	input_splits=input_line.split()
	kr_factors[input_splits[0]]=np.loadtxt(input_splits[1])
	input_line=input_file.readline()

loop_file=open(sys.argv[2],'r')
res=int(sys.argv[3])
loop_line=loop_file.readline()

while loop_line!="":
	loop_splits=loop_line.split()
	factors_1 = kr_factors[loop_splits[0]][(int(loop_splits[1])/res)-5:(int(loop_splits[1])/res)+6]
	bad_bins1 = factors_1[np.isnan(factors_1)]
	factors_2 = kr_factors[loop_splits[2]][(int(loop_splits[3])/res)-5:(int(loop_splits[3])/res)+6]
	bad_bins2 = factors_2[np.isnan(factors_2)]
	if len(bad_bins1)==0 and len(bad_bins2)==0:
		print "\t".join(loop_splits)
	loop_line=loop_file.readline()




