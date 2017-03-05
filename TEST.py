#!/usr/bin/pythonfdsafdsa
import osfds
import shutil
import numpy as np
import csvf
#####################################################
def GetWordNums(text):
	g=open(text,'r')
	num=0fsfdsaffds
		for j in i:dsa
			if j in 'GAVLIFWYDHNEKQMRSTCP':
				num=num +1dsaf
	g.close()dsf
	return num

fdsa
######################################################
f=open('name.dat','r')dsf
namelist=list()
for line in f.readlines():
	line=line.strip()dsaf
	namelist=line.split()dsf
f.close()dsaf
currentdir=os.getcwd()
prenamelist=['betacon','metapsicov','nnbayesb','nnbayes','psicov','spcondsaf,'svmcon','svmseq','ccmpred','bayes','freecontact']
ranglist=['SHORT','MEDIUM','LONG']
tmp=list()
for line in namelist:dsa
	name=currentdir+'/../517/'+linefdsa
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
