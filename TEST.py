#!/usr/bin/python
import os
import shutil
import numpy as np
import csv
#####################################################
def GetWordNums(text):
	g=open(text,'r')
	num=0
	for i in g.readlines():dsadsadsadsadd
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
ranglist=['SHORT','MEDIUM','LONG']
tmp=list()
for line in namelist:
	name=currentdir+'/../517/'+line
	os.chdir(name)
	length=GetWordNums('protein.seq')
	tmp.append(length)
lengthlist=np.array(tmp)
os.chdir(currentdir)



for i in range(517):
###############################################################################################
	mod=currentdir+'/../517/'+namelist[i]+'/native_contact_CA_all-reshape.dat'
	tmpmap=np.zeros((length,length))
	g=open(mod,'r')
	cc=0
	for k in g.readlines():
		cc+=1
	g.close()
	if(cc!=lengthlist[i]):
		print namelist[i]

##########################################true results###################################################################
