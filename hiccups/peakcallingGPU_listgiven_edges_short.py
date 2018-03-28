#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Very basic peak calling for GPU 
"""

import numpy as np
import sys
from scipy import sparse
from scipy import stats
from pycuda import driver, compiler, gpuarray, tools
import time
import copy
import array
import struct

# -- initialize the device
import pycuda.autoinit

kernel_code_template = """
__global__ void BasciPeakCallingKernel(float *c, float *expectedbl, float *expecteddonut, float *expectedh, float *expectedv, float *observed, float *b_bl, float *b_donut, float *b_h, float *b_v, float *p, float *tbl, float *td, float *th, float *tv, float *d, float *kr1, float *kr2, float *bound1, float *bound3)
{
    // 2D Thread ID 
    int t_col = threadIdx.x + blockIdx.x * blockDim.x;
    int t_row = threadIdx.y + blockIdx.y * blockDim.y;

    // Evalue is used to store the element of the matrix
    // that is computed by the thread


    float Evalue_bl =  0;
    float Edistvalue_bl = 0;
    float Evalue_donut =  0;
    float Edistvalue_donut = 0;
    float Evalue_h =  0;
    float Edistvalue_h = 0;
    float Evalue_v =  0;
    float Edistvalue_v = 0;
    float e_bl = 0;
    float e_donut = 0;
    float e_h = 0;
    float e_v = 0;
    float o = 0;
    float sbtrkt = 0;
    float bvalue_bl = 0;
    float bvalue_donut = 0;
    float bvalue_h = 0;
    float bvalue_v = 0;
    int wsize = %(WINDOW)s;
    int msize = %(MATRIX_SIZE)s;
    int pwidth = %(PEAKWIDTH)s;
    //int dvisor = %(DIVISOR)s;
    int diff = bound1[0] - bound3[0];
    while (abs(t_row+diff-t_col)<=(2*wsize)) {
	wsize = wsize - 1;
    }
    if (wsize<=pwidth) {
	wsize = pwidth + 1;
    }

	
    if (t_row>=20&&t_row<=(msize-20)&&t_col>=20&&t_col<=(msize-20)) {
    // calculate initial bottom left box
    for (int i = max(0,t_row+1); i < min(t_row+wsize+1, msize); ++i) {
	int test=max(0,t_col-wsize);
        for (int j = test; j < min(t_col, msize); ++j) {
      	    if (!isnan(c[i * msize + j])) { 	
		if (i+diff-j<0) {
		Evalue_bl += c[i * msize +j];
	  	Edistvalue_bl += d[abs(i+diff-j)];
                }
	    }
	}
    }
    //Subtract off the middle peak
    for (int i = max(0,t_row+1); i < min(t_row+pwidth+1, msize); ++i) {
	int test=max(0,t_col-pwidth);
        for (int j = test; j < min(t_col, msize); ++j) {
      	    if (!isnan(c[i * msize + j])) { 	
		if (i+diff-j<0) {
            	Evalue_bl -= c[i * msize +j];
	    	Edistvalue_bl -= d[abs(i+diff-j)];
                }
	    }
	}
    }
    
    //fix box dimensions
    while (Evalue_bl<16) {
	Evalue_bl=0;
	Edistvalue_bl=0;
	wsize+=1;
	//dvisor = powf(wsize,2.0) - powf(pwidth,2.0);
	for (int i = max(0,t_row+1); i < min(t_row+wsize+1, msize); ++i) {
		int test=max(0,t_col-wsize);
        	for (int j = test; j < min(t_col, msize); ++j) {
      	    	    if (!isnan(c[i * msize + j])) { 	
			if (i+diff-j<0) {
			Evalue_bl += c[i * msize +j];
	    	        Edistvalue_bl += d[abs(i+diff-j)];
        	    	if (i>=t_row+1) {
				if (i<t_row+pwidth+1) {
					if (j>=t_col-pwidth) {
						if (j<t_col) {
		    				Evalue_bl -= c[i * msize +j];
	    	    				Edistvalue_bl -= d[abs(i+diff-j)];
						}
					}
				}
			}
			}
		    }
		}
    	}

    	//Subtact off the middle peak
    	//for (int i = max(0,t_row+1); i < min(t_row+pwidth+1, msize); ++i) {
	//    int test=max(0,t_col-pwidth);
        //    for (int j = test; j < min(t_col, msize); ++j) {
      	//    	if (!isnan(c[i * msize + j])) { 	
		    //if (i+diff-j<0) {
	//	    Evalue_bl -= c[i * msize +j];
	//    	    Edistvalue_bl -= d[abs(i+diff-j)];
    	    	    //}
	//	}
	//    }
    	//}
	if (wsize == 20) {
		break;
	}
    	if (2*wsize>=abs(t_row+diff-t_col)) {
		break;
	}
    }
    
    // calculate donut
    for (int i = max(0,t_row-wsize); i < min(t_row+wsize+1, msize); ++i) {
	int test=max(0,t_col-wsize);
        for (int j = test; j < min(t_col+wsize+1, msize); ++j) {
      	    if (!isnan(c[i * msize + j])) { 	
		if (i+diff-j<0) {
      	  	Evalue_donut += c[i * msize +j];
	  	Edistvalue_donut += d[abs(i+diff-j)];
                }
	    }
	}
    }
    //Subtract off the middle peak
    for (int i = max(0,t_row-pwidth); i < min(t_row+pwidth+1, msize); ++i) {
	int test=max(0,t_col-pwidth);
        for (int j = test; j < min(t_col+pwidth+1, msize); ++j) {
      	    if (!isnan(c[i * msize + j])) { 	
		if (i+diff-j<0) {
            	Evalue_donut -= c[i * msize +j];
	  	Edistvalue_donut -= d[abs(i+diff-j)];
                }
	    }
	}
    }
    //Subtract off the cross hairs
    if ((t_row-pwidth)>0) {
    	for (int i = max(0,t_row-wsize); i < (t_row-pwidth); ++i) {
      	    if (!isnan(c[i * msize + t_col])) { 	
    		Evalue_donut -= c[i * msize + t_col];
		Edistvalue_donut -= d[abs(i+diff-t_col)];
	    }
	    for (int j = -1; j <=1 ; ++j) {
		Evalue_v += c[i * msize + t_col + j];
		Edistvalue_v += d[abs(i+diff-t_col-j)];
	    }
	}
    }
    if ((t_row+pwidth)<msize) {
    	for (int i = (t_row+pwidth+1); i < min(t_row+wsize+1,msize); ++i) {
      	    if (!isnan(c[i * msize + t_col])) { 	
    		Evalue_donut -= c[i * msize + t_col];
		Edistvalue_donut -= d[abs(i+diff-t_col)];
	    }
	    for (int j = -1; j <=1 ; ++j) {
		Evalue_v += c[i * msize + t_col + j];
		Edistvalue_v += d[abs(i+diff-t_col-j)];
	    }
	}
    }
    if ((t_col-pwidth)>0) {
    	for (int j = max(0,t_col-wsize); j < (t_col-pwidth); ++j) {
      	    if (!isnan(c[t_row * msize + j])) { 	
    		Evalue_donut -= c[t_row * msize + j];
		Edistvalue_donut -= d[abs(t_row+diff-j)];
	    }
	    for (int i = -1; i <=1 ; ++i) {
		Evalue_h += c[(t_row+i) * msize + j];
		Edistvalue_h += d[abs(t_row+i+diff-j)];
	    }
	}
    }
    if ((t_col+pwidth)<msize) {
    	for (int j = (t_col+pwidth+1); j < min(t_col+wsize+1,msize); ++j) {
      	    if (!isnan(c[t_row * msize + j])) { 	
    		Evalue_donut -= c[t_row * msize + j];
		Edistvalue_donut -= d[abs(t_row+diff-j)];
	    }
	    for (int i = -1; i <=1 ; ++i) {
		Evalue_h += c[(t_row+i) * msize + j];
		Edistvalue_h += d[abs(t_row+i+diff-j)];
	    }
	}
    }
    }
    //if (t_row+diff-t_col<(-1*pwidth)-2) {
    e_bl = ((Evalue_bl*d[abs(t_row+diff-t_col)])/Edistvalue_bl)*kr1[t_row]*kr2[t_col];
    e_donut = ((Evalue_donut*d[abs(t_row+diff-t_col)])/Edistvalue_donut)*kr1[t_row]*kr2[t_col];
    e_h = ((Evalue_h*d[abs(t_row+diff-t_col)])/Edistvalue_h)*kr1[t_row]*kr2[t_col];
    e_v = ((Evalue_v*d[abs(t_row+diff-t_col)])/Edistvalue_v)*kr1[t_row]*kr2[t_col];
    if (!isnan(e_bl)) {
	if (e_bl<=1) {
		bvalue_bl = 0;
	}
	else {
		bvalue_bl = floorf(logf(e_bl)/logf(powf(2.0,.33)));
	}
    }
    if (!isnan(e_donut)) {
	if (e_donut<=1) {
		bvalue_donut = 0;
	}
	else {
		bvalue_donut = floorf(logf(e_donut)/logf(powf(2.0,.33)));
	}
    }
    if (!isnan(e_h)) {
	if (e_h<=1) {
		bvalue_h = 0;
	}
	else {
		bvalue_h = floorf(logf(e_h)/logf(powf(2.0,.33)));
	}
    }
    if (!isnan(e_v)) {
	if (e_v<=1) {
		bvalue_v = 0;
	}
	else {
		bvalue_v = floorf(logf(e_v)/logf(powf(2.0,.33)));
	}
    }
    	
    // Write the matrix to device memory;
    // each thread writes one element
    expectedbl[t_row * msize + t_col] = e_bl;
    expecteddonut[t_row * msize + t_col] = e_donut;
    expectedh[t_row * msize + t_col] = e_h;
    expectedv[t_row * msize + t_col] = e_v;
    o = roundf(c[t_row * msize + t_col]*kr1[t_row]*kr2[t_col]);
    observed[t_row * msize + t_col] = o; //roundf(c[t_row * msize + t_col]*kr1[t_row]*kr2[t_col]);
    b_bl[t_row * msize + t_col] = bvalue_bl;
    b_donut[t_row * msize + t_col] = bvalue_donut;
    b_h[t_row * msize + t_col] = bvalue_h;
    b_v[t_row * msize + t_col] = bvalue_v;
    sbtrkt = fmaxf(tbl[(int) bvalue_bl],td[(int) bvalue_donut]);
    sbtrkt = fmaxf(sbtrkt, th[(int) bvalue_h]);
    sbtrkt = fmaxf(sbtrkt, tv[(int) bvalue_v]);
    p[t_row * msize + t_col] = o-sbtrkt;
}
"""
begin_time=time.time()

# take in input sheet
outfile1=open(sys.argv[2],'w')
outfile2=open(sys.argv[3],'w')
fdr=int(sys.argv[4])
targetlist=np.loadtxt(sys.argv[5])
hist_bl=np.zeros((40,10000))
hist_donut=np.zeros((40,10000))
hist_h=np.zeros((40,10000))
hist_v=np.zeros((40,10000))
fdrlog_bl=np.zeros((40,10000))
fdrlog_donut=np.zeros((40,10000))
fdrlog_h=np.zeros((40,10000))
fdrlog_v=np.zeros((40,10000))
threshold_bl=np.zeros(40, dtype=np.float32)
threshold_donut=np.zeros(40, dtype=np.float32)
threshold_h=np.zeros(40, dtype=np.float32)
threshold_v=np.zeros(40, dtype=np.float32)
bound1array=np.zeros(1, dtype=np.float32)
bound3array=np.zeros(1, dtype=np.float32)

MATRIX_SIZE = 540 #assume square for now
PEAKWIDTH = int(sys.argv[6])
WINDOW = int(sys.argv[7])
DIVISOR = (WINDOW)**2-(PEAKWIDTH)**2
#GPU stuff
#number of threads in block
block_size = 16
block =  block_size, block_size, 1 #threads per block = block_size*block_size
grid = int(MATRIX_SIZE/block_size)+1, int(MATRIX_SIZE/block_size)+1 # for grid of blocks
		
#print "Using block", block, " and grid", grid
		
# get the kernel code from the template 
# by specifying the constants
kernel_code = kernel_code_template % {
		  'MATRIX_SIZE': MATRIX_SIZE,
		  'PEAKWIDTH': PEAKWIDTH,
		  'WINDOW': WINDOW,
		  'DIVISOR': DIVISOR,
   	}
		
# compile the kernel code 
mod = compiler.SourceModule(kernel_code)
		
# get the kernel function from the compiled module
matrixmul = mod.get_function("BasciPeakCallingKernel")

for run in range(2):
	input_file=open(sys.argv[1],'r')
	input_line=input_file.readline()
	while input_line!="":
		input_splits=input_line.split()
		start_time=time.time()
	
		# load inputs
		f=open(input_splits[1],'rb')
		index1=array.array('i')
		index2=array.array('i')
		val=array.array('f')
		while True:
			try:
				index1.fromfile(f,1)
				index2.fromfile(f,1)
				val.fromfile(f,1)
			except EOFError: break
		val=map(np.float32,val)
		index1 = np.asarray(index1)
		index2 = np.asarray(index2)
		val = np.asarray(val)
		d_cpu=np.loadtxt(input_splits[2]).astype(np.float32)
		kr_total_cpu=np.loadtxt(input_splits[3]).astype(np.float32)
		res=int(input_splits[4])
		#threshold_bl_cpu=copy.deepcopy(threshold_bl)
		#threshold_donut_cpu=copy.deepcopy(threshold_donut)
		#print threshold_bl
		res=int(input_splits[4])
		
		load_time=time.time()
		print "Time to load chr"+input_splits[0]+" matrix: "+str(load_time-start_time)+"s"
		
		# make dense matrix on CPU memory
		#a_dense_cpu=np.asarray(sparse.coo_matrix((val,(index1,index2))).todense())
		#convert_time=time.time()
		#print "Time to convert chr"+input_splits[0]+" matrix to dense: "+str(convert_time-load_time)+"s"
		
		dim=np.amax(index1)
		dim2=np.amax(index2)
		for i in xrange((dim/500)+1):
			bound1_r=max((i*500),0)
			bound2_r=min((i+1)*500,dim)
			bound1=max((i*500)-20,0)
			bound2=min(((i+1)*500)+20,dim)
			if bound1==0:
				bound2=540
			if bound2==dim:
				bound1=dim-540
			diff1=bound1_r-bound1
			diff2=bound2-bound2_r
			cut1 = index1 >= bound1
			cut1_index1 = index1[cut1]
			cut1_index2 = index2[cut1]
			cut1_val = val[cut1]
			cut2 = cut1_index1 < bound2
			cut2_index1 = cut1_index1[cut2]
			cut2_index2 = cut1_index2[cut2]
			cut2_val = cut1_val[cut2] 
			for j in xrange(i,(dim/500)+1):
				bound3_r=max((j*500),0)
				bound4_r=min((j+1)*500,dim)
				bound3=max((j*500)-20,0)
				bound4=min(((j+1)*500)+20,dim)
				if bound3==0:
					bound4=540
				if bound4==dim:
					bound3=dim-540		
				diff3=bound3_r-bound3
				diff4=bound4-bound4_r
				cut3 = cut2_index2 >= bound3
				cut3_index1 = cut2_index1[cut3]
				cut3_index2 = cut2_index2[cut3]
				cut3_val = cut2_val[cut3]
				cut4 = cut3_index2 < bound4
				cut4_index1 = cut3_index1[cut4] - bound1
				cut4_index2 = cut3_index2[cut4] - bound3
				cut4_val = cut3_val[cut4]
				size_test1 = cut4_index1 == 539
				size_test2 = cut4_index2[size_test1] == 539
				if len(cut4_val[size_test2])==0:
					cut4_index1 = np.concatenate([cut4_index1,np.asarray([539])])
					cut4_index2 = np.concatenate([cut4_index2,np.asarray([539])])
					cut4_val = np.concatenate([cut4_val,np.asarray([0])])
				a_cpu=np.asarray(sparse.coo_matrix((cut4_val,(cut4_index1,cut4_index2))).todense(), np.float32)
				kr1_cpu=copy.deepcopy(kr_total_cpu[bound1:bound2])
				kr2_cpu=copy.deepcopy(kr_total_cpu[bound3:bound4])
				bound1array[0]=bound1	
				bound3array[0]=bound3
				
				gpu_time1=time.time()
				# transfer host (CPU) memory to device (GPU) memory 
				a_gpu = gpuarray.to_gpu(a_cpu) 
				d_gpu = gpuarray.to_gpu(d_cpu)
				kr1_gpu = gpuarray.to_gpu(kr1_cpu)
				kr2_gpu = gpuarray.to_gpu(kr2_cpu)
				threshold_bl_gpu = gpuarray.to_gpu(threshold_bl)
				threshold_donut_gpu = gpuarray.to_gpu(threshold_donut)
				threshold_h_gpu = gpuarray.to_gpu(threshold_h)
				threshold_v_gpu = gpuarray.to_gpu(threshold_v)
				bound1array_gpu = gpuarray.to_gpu(bound1array)
				bound3array_gpu = gpuarray.to_gpu(bound3array)
		
				# create empty gpu array for the result 
				expected_bl_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				expected_donut_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				expected_h_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				expected_v_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				observed_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				bin_bl_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				bin_donut_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				bin_h_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				bin_v_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
				peak_gpu = gpuarray.empty((np.shape(a_cpu)[0], np.shape(a_cpu)[1]), np.float32)
		
				# call the kernel on the card
				matrixmul(
					# inputs
					a_gpu, 
   				        # output
 				    	expected_bl_gpu,
					expected_donut_gpu,
 				    	expected_h_gpu,
					expected_v_gpu,
 				    	observed_gpu,
    				    	bin_bl_gpu,
					bin_donut_gpu,
    				    	bin_h_gpu,
					bin_v_gpu,
					peak_gpu,
					# thresholds
					threshold_bl_gpu,
					threshold_donut_gpu,
					threshold_h_gpu,
					threshold_v_gpu,
					# distance expected
   				    	d_gpu,
   		    			# kr
   		    			kr1_gpu,
   		    			kr2_gpu,
					# bounds
					bound1array_gpu,
					bound3array_gpu,
   		    			#  grid of blocks
    	            			grid = grid, 
    		    			#  block of threads
   		    			block = block
       		    		 	)
		
				expected_bl_result=expected_bl_gpu.get()
				expected_donut_result=expected_donut_gpu.get()
				expected_h_result=expected_h_gpu.get()
				expected_v_result=expected_v_gpu.get()
				observed_result=observed_gpu.get()
				bin_bl_result=bin_bl_gpu.get()
				bin_donut_result=bin_donut_gpu.get()
				bin_h_result=bin_h_gpu.get()
				bin_v_result=bin_v_gpu.get()
				peak_result=peak_gpu.get()
				gpu_time2=time.time()
				#print gpu_time2-gpu_time1, bound1_r, bound3_r
				observed_dense_cpu=np.zeros((540-diff2-diff1,540-diff4-diff3),dtype=int)
				observed_dense_cpu[:,:]=copy.deepcopy(observed_result[diff1:(540-diff2),diff3:(540-diff4)])
				#print observed_dense_cpu[observed_dense_cpu[:,:]>0]
				expected_bl_dense_cpu=copy.deepcopy(expected_bl_result[diff1:(540-diff2),diff3:(540-diff4)])
				expected_donut_dense_cpu=copy.deepcopy(expected_donut_result[diff1:(540-diff2),diff3:(540-diff4)])
				expected_h_dense_cpu=copy.deepcopy(expected_h_result[diff1:(540-diff2),diff3:(540-diff4)])
				expected_v_dense_cpu=copy.deepcopy(expected_v_result[diff1:(540-diff2),diff3:(540-diff4)])
				bin_bl_dense_cpu=copy.deepcopy(bin_bl_result[diff1:(540-diff2),diff3:(540-diff4)])
				bin_donut_dense_cpu=copy.deepcopy(bin_donut_result[diff1:(540-diff2),diff3:(540-diff4)])
				bin_h_dense_cpu=copy.deepcopy(bin_h_result[diff1:(540-diff2),diff3:(540-diff4)])
				bin_v_dense_cpu=copy.deepcopy(bin_v_result[diff1:(540-diff2),diff3:(540-diff4)])
				peak_dense_cpu=copy.deepcopy(peak_result[diff1:(540-diff2),diff3:(540-diff4)])
				if run==0:
					nanvals_bl = np.isnan(expected_bl_dense_cpu)
					nanvals_donut = np.isnan(expected_donut_dense_cpu)
					nanvals_h = np.isnan(expected_h_dense_cpu)
					nanvals_v = np.isnan(expected_v_dense_cpu)
					bin_bl_dense_cpu[nanvals_bl]=float('NaN')
					bin_bl_dense_cpu[nanvals_donut]=float('NaN')
					bin_bl_dense_cpu[nanvals_h]=float('NaN')
					bin_bl_dense_cpu[nanvals_v]=float('NaN')
					bin_donut_dense_cpu[nanvals_bl]=float('NaN')
					bin_donut_dense_cpu[nanvals_donut]=float('NaN')
					bin_donut_dense_cpu[nanvals_h]=float('NaN')
					bin_donut_dense_cpu[nanvals_v]=float('NaN')
					bin_h_dense_cpu[nanvals_bl]=float('NaN')
					bin_h_dense_cpu[nanvals_donut]=float('NaN')
					bin_h_dense_cpu[nanvals_h]=float('NaN')
					bin_h_dense_cpu[nanvals_v]=float('NaN')
					bin_v_dense_cpu[nanvals_bl]=float('NaN')
					bin_v_dense_cpu[nanvals_donut]=float('NaN')
					bin_v_dense_cpu[nanvals_h]=float('NaN')
					bin_v_dense_cpu[nanvals_v]=float('NaN')
					d_correct=(bound1_r-bound3_r)+int(sys.argv[6])+2
					dim_box=np.shape(observed_dense_cpu)[0]
					if d_correct>=(-1*dim_box):
						temp_cpu=np.ones(np.shape(observed_dense_cpu))
						nan_cpu=np.tril(temp_cpu,d_correct)
						nanvals_temp = nan_cpu==1
						nan_cpu[nanvals_temp]=float('NaN')
						del temp_cpu
						bin_bl_dense_cpu=bin_bl_dense_cpu+nan_cpu
						bin_donut_dense_cpu=bin_donut_dense_cpu+nan_cpu
						bin_h_dense_cpu=bin_h_dense_cpu+nan_cpu
						bin_v_dense_cpu=bin_v_dense_cpu+nan_cpu
					for k in xrange(40):
						blvals = bin_bl_dense_cpu == k
						donutvals = bin_donut_dense_cpu == k
						hvals = bin_h_dense_cpu == k
						vvals = bin_v_dense_cpu == k
						#hist, bin_edges = np.histogram(observed_dense_cpu[blvals],bins=range(10001))
						#hist_bl[i,:]=hist_bl[i,:]+hist
						hist_bl[k,:]=hist_bl[k,:]+np.bincount(observed_dense_cpu[blvals],minlength=10000)[0:10000]
						#hist, bin_edges = np.histogram(observed_dense_cpu[donutvals],bins=range(10001))
						hist_donut[k,:]=hist_donut[k,:]+np.bincount(observed_dense_cpu[donutvals],minlength=10000)[0:10000]
						#hist_donut[i,:]=hist_donut[i,:]+hist
						hist_h[k,:]=hist_h[k,:]+np.bincount(observed_dense_cpu[hvals],minlength=10000)[0:10000]
						hist_v[k,:]=hist_v[k,:]+np.bincount(observed_dense_cpu[vvals],minlength=10000)[0:10000]
				if run==1:
					if input_splits[0]=="X":
						chr=23
					else:
						chr=int(input_splits[0])
					test1=targetlist[:,0]==chr
					test2=targetlist[:,1]>=bound1_r*res
					test3=targetlist[:,1]<bound2_r*res
					test4=targetlist[:,4]>=bound3_r*res
					test5=targetlist[:,4]<bound4_r*res
					test6=targetlist[:,3]==chr
					inframe_peaks=targetlist[test1&test2&test3&test4&test5&test6]
					for k in xrange(len(inframe_peaks)):
						anchor1=inframe_peaks[k,1]/res-bound1_r
						anchor2=inframe_peaks[k,4]/res-bound3_r
						op=observed_dense_cpu[anchor1,anchor2]
						epbl=expected_bl_dense_cpu[anchor1,anchor2]
						epdonut=expected_donut_dense_cpu[anchor1,anchor2]
						eph=expected_h_dense_cpu[anchor1,anchor2]
						epv=expected_v_dense_cpu[anchor1,anchor2]
						bpbl=bin_bl_dense_cpu[anchor1,anchor2]
						bpdonut=bin_donut_dense_cpu[anchor1,anchor2]
						bph=bin_h_dense_cpu[anchor1,anchor2]
						bpv=bin_v_dense_cpu[anchor1,anchor2]
						if int(bpbl)<40 and op<10000 and int(bpdonut)<40 and int(bph)<40 and int(bpv)<40 and op>=0:
							fdr_bl=fdrlog_bl[int(bpbl),op]
							fdr_donut=fdrlog_donut[int(bpdonut),op]
							fdr_h=fdrlog_h[int(bph),op]
							fdr_v=fdrlog_v[int(bpv),op]
							outfile2.write(input_splits[0]+"\t"+str(inframe_peaks[k,1])+"\t"+input_splits[0]+"\t"+str(inframe_peaks[k,4])+"\t"+str(op)+"\t"+str(epbl)+"\t"+str(epdonut)+"\t"+str(eph)+"\t"+str(epv)+"\t"+str(bpbl)+"\t"+str(bpdonut)+"\t"+str(bph)+"\t"+str(bpv)+"\t"+str(fdr_bl)+"\t"+str(fdr_donut)+"\t"+str(fdr_h)+"\t"+str(fdr_v)+"\n")
		if run==0:
			hist_time=time.time()
			print "Time to calculate chr"+input_splits[0]+" expecteds and add to hist: "+str(hist_time-load_time)+"s"
		if run==1:	
			peak_time=time.time()
			print "Time to print chr"+input_splits[0]+" peaks: "+str(peak_time-load_time)+"s"	
		input_line=input_file.readline()
	input_file.close()
	if run==0:
		for i in xrange(40):#np.shape(observed_dense_cpu)[0]):
			expectedbl=sum(hist_bl[i,:])*stats.poisson.pmf(range(10000),2**((i+1)/3.0))
			expecteddonut=sum(hist_donut[i,:])*stats.poisson.pmf(range(10000),2**((i+1)/3.0))
			expectedh=sum(hist_h[i,:])*stats.poisson.pmf(range(10000),2**((i+1)/3.0))
			expectedv=sum(hist_v[i,:])*stats.poisson.pmf(range(10000),2**((i+1)/3.0))
			if sum(hist_bl[i,:])>0:
				for j in xrange(10000):
					if fdr*sum(expectedbl[j:10000])<=sum(hist_bl[i,j:10000]):
						threshold_bl[i]=np.float32(j-1)
						break
				for j in xrange(10000):
					sum1=sum(expectedbl[j:10000])
					sum2=sum(hist_bl[i,j:10000])
					if sum2>0:
						fdrlog_bl[i,j]=sum1/(sum2*1.0)
					else:
						break
			else:
				threshold_bl[i]=np.float32(10000)
			if sum(hist_donut[i,:])>0:
				for j in xrange(10000):
					if fdr*sum(expecteddonut[j:10000])<=sum(hist_donut[i,j:10000]):
						threshold_donut[i]=np.float32(j-1)
						break
				for j in xrange(10000):
					sum1=sum(expecteddonut[j:10000])
					sum2=sum(hist_donut[i,j:10000])
					if sum2>0:
						fdrlog_donut[i,j]=sum1/(sum2*1.0)
					else:
						break
			else:
				threshold_donut[i]=np.float32(10000)
			if sum(hist_h[i,:])>0:
				for j in xrange(10000):
					if fdr*sum(expectedh[j:10000])<=sum(hist_h[i,j:10000]):
						threshold_h[i]=np.float32(j-1)
						break
				for j in xrange(10000):
					sum1=sum(expectedh[j:10000])
					sum2=sum(hist_h[i,j:10000])
					if sum2>0:
						fdrlog_h[i,j]=sum1/(sum2*1.0)
					else:
						break
			else:
				threshold_h[i]=np.float32(10000)
			if sum(hist_v[i,:])>0:
				for j in xrange(10000):
					if fdr*sum(expectedv[j:10000])<=sum(hist_v[i,j:10000]):
						threshold_v[i]=np.float32(j-1)
						break
				for j in xrange(10000):
					sum1=sum(expectedv[j:10000])
					sum2=sum(hist_v[i,j:10000])
					if sum2>0:
						fdrlog_v[i,j]=sum1/(sum2*1.0)
					else:
						break
			else:
				threshold_v[i]=np.float32(10000)
			outfile1.write(str(i)+"\t"+str(threshold_bl[i])+"\t"+str(threshold_donut[i])+"\t"+str(threshold_h[i])+"\t"+str(threshold_v[i])+"\n")		
		thresh_time=time.time()
		print "Time to calculate thresholds: "+str(thresh_time-hist_time)+"s"	
final_time=time.time()
print "Total time: "+str(final_time-begin_time)+"s"	
outfile1.close()
outfile2.close()

