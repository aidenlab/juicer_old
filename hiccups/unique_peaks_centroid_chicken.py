#!/usr/bin/env python

from numpy import *
import sys

file=sys.argv[1]
d=int(sys.argv[2])
initial_d=int(sys.argv[2])
res=int(sys.argv[3])
input_color=sys.argv[4]
fdr_color_switch=sys.argv[5]
peaklist=loadtxt(file)
newpeaklist=zeros(shape=(1,17))
r=0
while len(peaklist[:,0])!=0:
	#print amin(peaklist[:,7])
	peakindex=argmax(peaklist[:,4])
	#print peaklist[peakindex,:]
	chr=peaklist[peakindex,0]
	xcoord=peaklist[peakindex,1]
	ycoord=peaklist[peakindex,3]
	newpeaklist=peaklist[peakindex,:]
	deletelist=[]
	deletelist.append(peakindex)
	xcoordlist=[]
	ycoordlist=[]
	xcoordlist.append(peaklist[peakindex,1])
	ycoordlist.append(peaklist[peakindex,3])
	#print len(peaklist[:,0])
	for i in range(len(peaklist[:,0])):
		if peaklist[i,0]==chr and sqrt((peaklist[i,1]-xcoord)**2+(peaklist[i,3]-ycoord)**2)<=d and i!=peakindex:
			#print peaklist[peakindex,:], peaklist[i,:]
			#print i, peaklist[i,:]
			deletelist.append(i)
			xcoordlist.append(peaklist[i,1])
			ycoordlist.append(peaklist[i,3])
			xcoord=mean(xcoordlist)
			ycoord=mean(ycoordlist)
			r=0
			for j in range(len(xcoordlist)):
				rprime=sqrt((xcoordlist[j]-xcoord)**2+(ycoordlist[j]-ycoord)**2)
				#print i,rprime
				if rprime>r:
					r=rprime
			d=initial_d+r
	#newpeaklist[1]=xcoord
	#newpeaklist[3]=ycoord
	deletelist.sort(reverse=True)
	for j in deletelist:		
		#print j
		peaklist=delete(peaklist,(j),axis=0)
	numcollapsed=len(deletelist)

	color=input_color
	if fdr_color_switch=='y':
		fdr=-1.0*floor(log10(max(newpeaklist[9],newpeaklist[10])))
		if fdr>10:
			fdr=10
		color_splits=input_color.split(",")
		red=str(int(int(color_splits[0])*(fdr/10)))
		green=str(int(int(color_splits[1])*(fdr/10)))
		blue=str(int(int(color_splits[2])*(fdr/10)))
		color=red+","+green+","+blue 
	if newpeaklist[0]==29:
		print "chrZ"+"\t"+str(int(newpeaklist[1]))+"\t"+str(int(newpeaklist[1])+res)+"\t"+"chrZ"+"\t"+str(int(newpeaklist[3]))+"\t"+str(int(newpeaklist[3])+res)+"\t"+color+"\t"+str(int(newpeaklist[4]))+"\t"+str(newpeaklist[5])+"\t"+str(newpeaklist[6])+"\t"+str(newpeaklist[7])+"\t"+str(newpeaklist[8])+"\t"+str(int(newpeaklist[9]))+"\t"+str(int(newpeaklist[10]))+"\t"+str(int(newpeaklist[11]))+"\t"+str(int(newpeaklist[12]))+"\t"+str(newpeaklist[13])+"\t"+str(newpeaklist[14])+"\t"+str(newpeaklist[15])+"\t"+str(newpeaklist[16])+"\t"+str(int(numcollapsed))+"\t"+str(int(xcoord+(res/2)))+"\t"+str(int(ycoord+(res/2)))+"\t"+str(r)
	else:
		print "chr"+str(int(newpeaklist[0]))+"\t"+str(int(newpeaklist[1]))+"\t"+str(int(newpeaklist[1])+res)+"\t"+"chr"+str(int(newpeaklist[2]))+"\t"+str(int(newpeaklist[3]))+"\t"+str(int(newpeaklist[3])+res)+"\t"+color+"\t"+str(int(newpeaklist[4]))+"\t"+str(newpeaklist[5])+"\t"+str(newpeaklist[6])+"\t"+str(newpeaklist[7])+"\t"+str(newpeaklist[8])+"\t"+str(int(newpeaklist[9]))+"\t"+str(int(newpeaklist[10]))+"\t"+str(int(newpeaklist[11]))+"\t"+str(int(newpeaklist[12]))+"\t"+str(newpeaklist[13])+"\t"+str(newpeaklist[14])+"\t"+str(newpeaklist[15])+"\t"+str(newpeaklist[16])+"\t"+str(int(numcollapsed))+"\t"+str(int(xcoord+(res/2)))+"\t"+str(int(ycoord+(res/2)))+"\t"+str(r)
	r=0
