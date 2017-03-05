#!/usr/bin/python
import os
import shutil
import numpy as np
import csv
#####################################################
def GetWordNums(text):
	g=open(text,'r')
	num=0
	for i in g.readlines():
		for j in i:
			if j in 'GAVLIFWYDHNEKQMRSTCP':
				num=num +1
	g.close()
	return num


######################################################
f=open('name.dat','r')
namelist=list()
for line in f.readlines():
	line=line.strip()
	namelist=line.split()
f.close()
currentdir=os.getcwd()
prenamelist=['betacon','metapsicov','nnbayesb','nnbayes','psicov','spcon','svmcon','svmseq','ccmpred','bayes','freecontact']
ranglist=['SHORT']
tmp=list()
for line in namelist:
	name=currentdir+'/../517/'+line
	os.chdir(name)
	length=GetWordNums('protein.seq')
	tmp.append(length)
lengthlist=np.array(tmp)
os.chdir(currentdir)
#################################################################################################
def return_batch(NUM,ranges,windowsize,filename):#NUM:number of batch,range:SHORT,MEDIUM and LONG, windowsize:11 and so on. length:protein length,filnema=file name
	downlimit=6
	uplimit=10000
	NUM=int(NUM)
	windowsize=int(windowsize)
	if(ranges=='LONG'):
		downlimit=25
		uplimit=10000
	elif(ranges=='MEDIUM'):
		downlimit=12
		uplimit=24
	elif(ranges=='SHORT'):
		downlimit=6
		uplimit=11
	else:
		print 'ALL RANGE!'
##############################################################
	batch=np.zeros((windowsize,windowsize,11,NUM))
	native=np.zeros((NUM,2))
	SS=np.zeros((NUM,2*windowsize,3))
	SOLV=np.zeros((NUM,2*windowsize))
	COLSTATS=np.zeros((NUM,2*windowsize,21))
###############################################################
	for count in range(NUM):
		if(count%100==0):
			print count
##########################################contact map input#####################################################################
		i=0
		j=0
		proteinNo=np.random.randint(0,502)#the last 15 proteins are not used
		length=lengthlist[proteinNo]
		while(abs(i-j)<downlimit or abs(i-j)>uplimit):
			i=np.random.randint(0,length)
			j=np.random.randint(0,length)
		namecount=0
#		print namelist[proteinNo],length
			
		for name in prenamelist:
			mod=currentdir+'/../517/'+namelist[proteinNo]+'/'+name+'-reshape.dat'
			tmpmap=np.zeros((length,length))
			g=open(mod,'r')
			cc=0
			for k in g.readlines():
				tmpmap[cc,:]=np.array(map(eval,k.strip().split(',')),dtype='float64')
				cc+=1
			g.close()
			for k in range(windowsize):
				for l in range(windowsize):
					if i-int(windowsize/2)+k>=0 and i-int(windowsize/2)+k<length and j-int(windowsize/2)+l>=0 and j+int(windowsize/2)+l<length :
						batch[k,l,namecount,count]=tmpmap[i-int(windowsize)+k,j-int(windowsize)+l]
					else:
						batch[k,l,namecount,count]=0.0
##########################################true results###################################################################
		mod=currentdir+'/../517/'+namelist[proteinNo]+'/native_contact_CA_all-reshape.dat'
		tmpmap=np.zeros((length,length))
		g=open(mod,'r')
		cc=0
		for k in g.readlines():
			tmpmap[cc,:]=np.array(map(eval,k.strip().split(',')),dtype='float64')
			cc+=1
		g.close()
		native[count,0]=tmpmap[i,j]
		mod=currentdir+'/../517/'+namelist[proteinNo]+'/native_contact_CB_all-reshape.dat'
		tmpmap=np.zeros((length,length))
		g=open(mod,'r')
		cc=0
		for k in g.readlines():
			tmpmap[cc,:]=np.array(map(eval,k.strip().split(',')),dtype='float64')
			cc+=1
		g.close()
		native[count,1]=tmpmap[i,j]
#####################################################Secondary structure############################################
		tmpmap=np.zeros((length,3))
		mod=currentdir+'/../517/'+namelist[proteinNo]+'/protein.ss'
		g=open(mod,'r')
		cc=0
		for k in g.readlines():
			tmpmap[cc,0]=float(k.strip().split()[3])
			tmpmap[cc,1]=float(k.strip().split()[4])
			tmpmap[cc,2]=float(k.strip().split()[5])
			cc+=1
		g.close()
		for k in range(windowsize):
			if i-int(windowsize/2)+k>=0 and i-int(windowsize/2)+k<length:
				SS[count,k,0]=tmpmap[i-int(windowsize/2)+k,0]
				SS[count,k,1]=tmpmap[i-int(windowsize/2)+k,1]
				SS[count,k,2]=tmpmap[i-int(windowsize/2)+k,2]
			else:
				SS[count,k,0]=0.0
				SS[count,k,1]=0.0
				SS[count,k,2]=0.0
		for k in range(windowsize):
			if j-int(windowsize/2)+k>=0 and j-int(windowsize/2)+k<length:
				SS[count,windowsize+k,0]=tmpmap[j-int(windowsize/2)+k,0]
				SS[count,windowsize+k,1]=tmpmap[j-int(windowsize/2)+k,1]
				SS[count,windowsize+k,2]=tmpmap[j-int(windowsize/2)+k,2]
			else:
				SS[count,windowsize+k,0]=0.0
				SS[count,windowsize+k,1]=0.0
				SS[count,windowsize+k,2]=0.0

		

###################################################solvent accessibility###########################################
		tmpmap=np.zeros(length)
		mod=currentdir+'/../517/'+namelist[proteinNo]+'/protein.solv'
		g=open(mod,'r')
		cc=0
		for k in g.readlines():
			tmpmap[cc]=float(k.strip().split()[2])
			cc+=1
		g.close()
		for k in range(windowsize):
			if i-int(windowsize/2)+k>=0 and i-int(windowsize/2)+k<length:
				SOLV[count,k]=tmpmap[i-int(windowsize/2)+k]
			else:
				SOLV[count,k]=0.0
		for k in range(windowsize):
			if j-int(windowsize/2)+k>=0 and j-int(windowsize/2)+k<length:
				SOLV[count,windowsize+k]=tmpmap[j-int(windowsize/2)+k]
			else:
				SOLV[count,windowsize+k]=0.0

###################################################colstats########################################################
		tmpmap=np.zeros((length,21))
		mod=currentdir+'/../517/'+namelist[proteinNo]+'/protein.colstats'
		g=open(mod,'r')
		g.readline()
		g.readline()
		g.readline()
		g.readline()
		cc=0
		for k in g.readlines():
			for l in range(21):
				tmpmap[cc][l]=float(k.strip().split()[l])
			cc+=1
		g.close()
		for k in range(windowsize):
			for m in range(21):
				if i-int(windowsize/2)+k>=0 and i-int(windowsize/2)+k<length:
					COLSTATS[count,k,m]=tmpmap[i-int(windowsize/2)+k,m]
				else:
					COLSTATS[count,k,m]=0.0
		for k in range(windowsize):
			for m in range(21):
				if j-int(windowsize/2)+k>=0 and j-int(windowsize/2)+k<length:
					COLSTATS[count,windowsize+k,m]=tmpmap[j-int(windowsize/2)+k,m]
				else:
					COLSTATS[count,windowsize+k,m]=0.0		
##########################################save files###############################################################
	filenames=currentdir+'/'+ranges+'/'+filename
	np.save(filenames,batch)
	filenames=currentdir+'/'+ranges+'/'+filename+'-Truelabel'
	np.save(filenames,native)
	filenames=currentdir+'/'+ranges+'/'+filename+'-SS'
	np.save(filenames,SS)
	filenames=currentdir+'/'+ranges+'/'+filename+'-SOLV'
	np.save(filenames,SOLV)
	filenames=currentdir+'/'+ranges+'/'+filename+'-COLSTATS'
	np.save(filenames,COLSTATS)
	
################################################################################################
for ranges in ranglist:
	for i in range(1000):
		print ranges
		filename=ranges+str(i)
		return_batch(10000,ranges,11,filename)
