#!/usr/bin/env python
# coding: utf-8
# In[ ]:
import sys
from math import *
import copy
import time
import random
n_cond = 23814 #number of contidions
n_moas = 206 #number of MoAs
n_gen = 772 #number of genes
n_cell = 100 #number of cell types
#lambdag=0.25 #contribution to the likelihood genes (aprox. #genes = 4*#MoAs) --> weight: 1/4 
#lambda g, Toñi defines it as 0.3*18.384408 ? Ask her!!
#lambdav=0.25 #contribution to the likelihood viability
num=sys.argv[1]	#when distributing the code in the cluster, this will be the ID of this random realization: integer number 0,1, etc
K=int(sys.argv[2]) #groups of conditions
L=int(sys.argv[3]) #groups of moas
T=int(sys.argv[4]) #groups of genes
runs=int(sys.argv[5]) #how many times the iterative algorithm will run
fold=int(sys.argv[6]) #once the dataset has been partitioned, we can use different folds (only affecting the moa-cond files)
lambdag=float(sys.argv[7])
direct="../Data/"
#read the data in the most efficient way: for the c-MoAs separate the data between 0 and 1
#there are three type of interaccions (with MoAs, genes and cells). Let's start by the first two: c-MoAs and c-genes
################################
##### Mechanisms of Action #####
################################
### Cond-MoA with 0s ###
fh2=open(direct+'c-moastrain0.csv','r')
igot2=fh2.readlines() #reading each row as a string, list of strings
cmoas0=[]
for line in igot2:
	about = line.strip().split(',') #list of strings (for each row)
	cmoas0.append((int(about[0]),int(about[1]))) #1st element: condition integer; 2nd element: MoA integer (for those combinations with a 0)
fh2.close()
### Cond-MoA with 1s ###
fh2=open(direct+'c-moastrain1.csv','r')
igot2=fh2.readlines()
cmoas1=[]
for line in igot2:
	about = line.strip().split(',')
	cmoas1.append((int(about[0]),int(about[1])))
fh2.close()
################################
########### Genes ##############
################################
### Cond-Gene with -1s ###
fh2=open(direct+'/c-genestrain-1.csv','r')
igot2=fh2.readlines()
cgenes_1=[]
for line in igot2:
	about = line.strip().split(',')
	cgenes_1.append((int(about[0]),int(about[1])))
fh2.close()
### Cond-Gene with 0s ###
fh2=open(direct+'c-genestrain0.csv','r')
igot2=fh2.readlines()
cgenes0=[]
for line in igot2:
	about = line.strip().split(',')
	cgenes0.append((int(about[0]),int(about[1])))
fh2.close()
### Cond-Gene with 1s ###
fh2=open(direct+'c-genestrain1.csv','r')
igot2=fh2.readlines()
cgenes1=[]
for line in igot2:
	about = line.strip().split(',')
	cgenes1.append((int(about[0]),int(about[1])))
fh2.close()
########################
##### Random param.#####
########################
## Theta:
theta=[]
ntheta=[]
dtheta=[]
for c in range(n_cond): 
	vec=[]
	for alpha in range(K): #K = number of conditions groups
		vec.append(random.random())
	theta.append(vec)
	ntheta.append([0.]*K) #same dimensions as theta but with 0.0 everywhere
	dtheta.append([0.]*K) #same dimensions as theta but with 0.0 everywhere
## Eta:
eta=[]
neta=[]
deta=[]
for m in range(n_moas):
	vec=[]
	for beta in range(L):
		vec.append(random.random())
	eta.append(vec)
	neta.append([0.]*L)
	deta.append([0.]*L)
## Phi:
phi=[]
nphi=[]
dphi=[]
for g in range(n_gen):
	vec=[]
	for gamma in range(T):
		vec.append(random.random())
	phi.append(vec)
	nphi.append([0.]*T)
	dphi.append([0.]*T)	
## P matrix:
pr=[]
npr=[]
dpr=[]
for alpha in range(K):
	vec=[]
	for beta in range(L):
		vec.append(random.random())
	pr.append(vec)
	npr.append([0.]*L)
	dpr.append([0.]*L)
## Q matrix: 3-dimensional (since we have 3 different link values, we can not do the same as with pmat) 
qr=[]
nqr=[]
dqr=[]
for e in range(3): #e can be -1, 0 or 1 --> in here, translated as index 0 (for link -1), 1 (for link 0), 2 (for link 1)
	qr.append([])
	nqr.append([])
	for alpha in range(K):
		qr[e].append([])
		nqr[e].append([])
		vec=[]
		for gamma in range(T): 
			vec.append(random.random())
		qr[e][alpha] = vec 
		nqr[e][alpha] = [0.]*T
for alpha in range(K):
    dqr.append([0.]*T) #dqr has not the same shape because it's where we'll save the sum of psi[alpha][beta] 
#for any e link --> for the normalization when computing the iterative algorithm
########################
#####Normalizations#####
########################
## Theta norm:
for c in range(n_cond): #each condition
	D=0.
	for alpha in range(K): #each condition group
		D=D+theta[c][alpha] #for each condition, D will be the sum of the elements of the membership vector
	for alpha in range(K):
		theta[c][alpha]=theta[c][alpha]/D #for each condition, divide all its elements by D
## Eta norm:
for m in range(n_moas):
	D=0.
	for beta in range(L):
		D=D+eta[m][beta]
	for beta in range(L):
		eta[m][beta]=eta[m][beta]/D
## Phi norm:
for g in range(n_gen):
	D=0.
	for gamma in range(T):
		D=D+phi[g][gamma]
	for gamma in range(T):
		phi[g][gamma]=phi[g][gamma]/D
## P norm: there's no point in doing it (since we don't use the dimension for specifying the value of the link,
## in this matrix we hace the probabilities(1), if we want the proba(0) we'll say that it is 1-proba(1). So,
## we already have the normalization impicitly (we're forcing that sum to be 1))
###################
## Q norm: this matrix contains the dimension mentioned above --> we must normalize it in order to have a result of 1 
## when adding all the probabilities for a fixed pair of cond-gene (but modifying the link value) to be 1.
for alpha in range(K): #condition groups
	for gamma in range(T): #gene groups
		D = 0.
		for e in range(3): #link
			D += qr[e][alpha][gamma]     
		for e in range(3):
			qr[e][alpha][gamma] = qr[e][alpha][gamma]/D
#############################
#############################
#####Iterative algorithm#####
#############################
#############################
Runs=int(sys.argv[5]) #total number of iterations for the local maximum likelihood. AG to Gemma: We should discuss this bit
Old=-1000000000000
Like=-1000000000
iteration=0
fout=open('likelihood_K'+str(K)+'L'+str(L)+'T'+str(T)+'_fold'+str(fold)+'_iters'+str(runs)+'_run'+str(num)+'.csv','w')
print('starting iterative algorithm...')
while (Old-Like)<0.1 and iteration<Runs:	#this should not happen, but we kill the process if the likelihood decrease
	print('iteration', iteration)
	start=time.time()
	Old=Like
	iteration=iteration+1
	#first we feed the model parameters with Cond-MoAs 1s
	for pair in cmoas1:
		condition=pair[0] 
		moa=pair[1] 
		D_omega=0.	#this is the denominator of the probability function omega in eq. 9
		for alpha in range(K):
			for beta in range(L):
				D_omega=D_omega+theta[condition][alpha]*eta[moa][beta]*pr[alpha][beta]
		for alpha in range(K):
			for beta in range(L):
				omega=theta[condition][alpha]*eta[moa][beta]*pr[alpha][beta]/D_omega	#auxiliary probability function eq.9 (omega)
				#(alpha*AA+CC)*
				ntheta[condition][alpha]=ntheta[condition][alpha]+omega #es van sumant els valors de les omegues (faltarà que ho dividim per Mc)
				neta[moa][beta]=neta[moa][beta]+omega
				npr[alpha][beta]=npr[alpha][beta]+omega #new p matrix (withou norm.), only interested in link 1
				dpr[alpha][beta]=dpr[alpha][beta]+omega #dpr will save all the values for the future normalization
	#Cond-MoAs with 0s
	for pair in cmoas0:
		condition=pair[0]
		moa=pair[1]
		D_omega=0.
		for alpha in range(K):
			for beta in range(L):
				D_omega=D_omega+theta[condition][alpha]*eta[moa][beta]*(1.-pr[alpha][beta])
		
		for alpha in range(K):
			for beta in range(L):
				omega=theta[condition][alpha]*eta[moa][beta]*(1.-pr[alpha][beta])/D_omega
				#(alpha*AA+CC)*
				ntheta[condition][alpha]=ntheta[condition][alpha]+omega
				neta[moa][beta]=neta[moa][beta]+omega
				dpr[alpha][beta]=dpr[alpha][beta]+omega #dpr will save all the values for the future normalization
	#Cond-Genes with -1s
	for pair in cgenes_1:
		condition=pair[0]
		gene=pair[1]
		D_psi=0.
		for alpha in range(K):
			for gamma in range(T):
				D_psi=D_psi+theta[condition][alpha]*phi[gene][gamma]*qr[0][alpha][gamma] #0 corresponds to link -1
		
		for alpha in range(K):
			for gamma in range(T):
				psi=theta[condition][alpha]*phi[gene][gamma]*qr[0][alpha][gamma]/D_psi
				#(alpha*AA+CC)*
				ntheta[condition][alpha]=ntheta[condition][alpha]+lambdag*psi
				nphi[gene][gamma]=nphi[gene][gamma]+psi
				nqr[0][alpha][gamma]=nqr[0][alpha][gamma]+psi             
				dqr[alpha][gamma]=dqr[alpha][gamma]+psi #for future normalization
	#Cond-Genes with 0s
	for pair in cgenes0:
		condition=pair[0]
		gene=pair[1]
		D_psi=0.
		for alpha in range(K):
			for gamma in range(T):
				D_psi=D_psi+theta[condition][alpha]*phi[gene][gamma]*qr[1][alpha][gamma] #1 corresponds to link 0
		
		for alpha in range(K):
			for gamma in range(T):
				psi=theta[condition][alpha]*phi[gene][gamma]*qr[1][alpha][gamma]/D_psi
				#(alpha*AA+CC)*
				ntheta[condition][alpha]=ntheta[condition][alpha]+lambdag*psi
				nphi[gene][gamma]=nphi[gene][gamma]+psi
				nqr[1][alpha][gamma]=nqr[1][alpha][gamma]+psi 
				dqr[alpha][gamma]=dqr[alpha][gamma]+psi 
	#Cond-Genes with 1s
	for pair in cgenes1:
		condition=pair[0]
		gene=pair[1]
		D_psi=0.
		for alpha in range(K):
			for gamma in range(T):
				D_psi=D_psi+theta[condition][alpha]*phi[gene][gamma]*qr[2][alpha][gamma] #2 corresponds to link 1
		
		for alpha in range(K):
			for gamma in range(T):
				psi=theta[condition][alpha]*phi[gene][gamma]*qr[2][alpha][gamma]/D_psi
				#(alpha*AA+CC)*
				ntheta[condition][alpha]=ntheta[condition][alpha]+lambdag*psi
				nphi[gene][gamma]=nphi[gene][gamma]+psi
				nqr[2][alpha][gamma]=nqr[2][alpha][gamma]+psi 
				dqr[alpha][gamma]=dqr[alpha][gamma]+psi 
	#normalization:
	## theta:  
	for condition in range(n_cond): 
		sumi=sum(ntheta[condition])
		for alpha in range(K):
			try:
				ntheta[condition][alpha]=ntheta[condition][alpha]/(sumi)
			except ZeroDivisionError:
				ntheta[condition][alpha]=0.
	## eta:      
	for moa in range(n_moas):
		sumi=sum(neta[moa])
		for beta in range(L):
			try:
				neta[moa][beta]=neta[moa][beta]/(sumi)
			except ZeroDivisionError:
				neta[moa][beta]=0.
	## phi:      
	for gene in range(n_gen):
		sumi=sum(nphi[gene])
		for gamma in range(T):
			try:
				nphi[gene][gamma]=nphi[gene][gamma]/(sumi)
			except ZeroDivisionError:
				nphi[gene][gamma]=0.
	## p matrix:                
	for alpha in range(K):
		for beta in range(L):
			npr[alpha][beta]=npr[alpha][beta]/(dpr[alpha][beta])
	## q matrix:     
	for alpha in range(K): #condition groups
		for gamma in range(T): #gene groups
			for e in range(3):
				nqr[e][alpha][gamma] = nqr[e][alpha][gamma]/(dqr[alpha][gamma])
		               
	#update new values of the parameters
	theta=copy.copy(ntheta)
	eta=copy.copy(neta)	
	phi=copy.copy(nphi)   
	for alpha in range(K):
		for beta in range(L):
			pr[alpha][beta]=npr[alpha][beta]
	for e in range(3):
		for alpha in range(K):
			for gamma in range(T):
				qr[e][alpha][gamma]=nqr[e][alpha][gamma]
	for condition in range(n_cond):
		ntheta[condition]=[0.]*K
	for moa in range(n_moas):
		neta[moa]=[0.]*L
	for gene in range(n_gen):
		nphi[gene]=[0.]*T       
	for alpha in range(K):
		for beta in range(L):
			npr[alpha][beta]=0.
			dpr[alpha][beta]=0.
	for alpha in range(K):
		for gamma in range(T):
			dqr[alpha][gamma]=0. #without link (normalization)
			for e in range(3):
				nqr[e][alpha][gamma]=0.
		            
	#compute the log posterior:
	log_posterior_moa = 0.
	log_posterior_gene = 0.
	#Cond-MoAs with 1s
	for pair in cmoas1:
		condition=pair[0]
		moa=pair[1]
		adding_moa=0.
		for alpha in range(K):
			for beta in range(L):
				proba=theta[condition][alpha]*eta[moa][beta]*pr[alpha][beta]
				adding_moa+=proba
		log_posterior_moa += log(adding_moa)
	#Cond-MoAs with 0s
	for pair in cmoas0:
		condition=pair[0]
		moa=pair[1]
		adding_moa=0.
		for alpha in range(K):
			for beta in range(L):
				proba=theta[condition][alpha]*eta[moa][beta]*(1.-pr[alpha][beta])
				adding_moa+=proba
		log_posterior_moa += log(adding_moa)
	#Cond-Gene with -1s
	for pair in cgenes_1:
		condition=pair[0]
		gene=pair[1]
		adding_gene=0.
		for alpha in range(K):
			for gamma in range(T):
				proba=theta[condition][alpha]*phi[gene][gamma]*qr[0][alpha][gamma]
				adding_gene+=proba
		log_posterior_gene += log(adding_gene)       
	#Cond-Gene with 0s
	for pair in cgenes0:
		condition=pair[0]
		gene=pair[1]
		adding_gene=0.
		for alpha in range(K):
			for gamma in range(T):
				proba=theta[condition][alpha]*phi[gene][gamma]*qr[1][alpha][gamma]
				adding_gene+=proba
		log_posterior_gene += log(adding_gene)
	#Cond-Gene with 1s
	for pair in cgenes1:
		condition=pair[0]
		gene=pair[1]
		adding_gene=0.
		for alpha in range(K):
			for gamma in range(T):
				proba=theta[condition][alpha]*phi[gene][gamma]*qr[2][alpha][gamma]
				adding_gene+=proba
		log_posterior_gene += log(adding_gene)
	Like = log_posterior_moa + log_posterior_gene*lambdag
	print('genes', log_posterior_gene, 'moa', log_posterior_moa)
	print(Like, time.time()-start)
	fout.write('%s,%f\n' % (iteration,Like))
fout.close()
#save parameters
fout2=open('param_K'+str(K)+'L'+str(L)+'T'+str(T)+'_fold'+str(fold)+'_iters'+str(runs)+'_run'+str(num)+'.dat','w') #num: integer (see first lines of code)
fout2.write('%s %s\n' % (iteration,Like))
for condition in range(n_cond):
	for alpha in range(K):
		fout2.write('%s ' % theta[condition][alpha])
	fout2.write('\n')
					
for moa in range(n_moas):
	for beta in range(L):
		fout2.write('%s ' % eta[moa][beta])
	fout2.write('\n')
					
for gene in range(n_gen):
	for gamma in range(T):
		fout2.write('%s ' % phi[gene][gamma])
	fout2.write('\n')
					
for alpha in range(K):
	for beta in range(L):
		fout2.write('%s ' % pr[alpha][beta])
	fout2.write('\n')
					
for e in range(3):
	for alpha in range(K):
		for gamma in range(T):
			fout2.write('%s ' % qr[e][alpha][gamma])
		fout2.write('\n')
	fout2.write('\n')
	
fout2.close()					
#Checking normalizations:
print('\nChecking normalizations:')
print('any wrong theta?')
for c in range(n_cond):
	sumi=sum(theta[c])
	if sumi<0.99 or sumi>1.01:
		print(c, ':', 'norm =', sumi)
print('any wrong eta?')
for m in range(n_moas):
	sumi=sum(eta[m])
	if sumi<0.99 or sumi>1.01:
		print(m, ':', 'norm =', sumi)
print('any wrong phi?')
for g in range(n_gen):
	sumi=sum(phi[g])
	if sumi<0.99 or sumi>1.01:
		print(g, ':', 'norm =', sumi)
print('any wrong qr?')
for alpha in range(K):
	for gamma in range(T):
		sumi=qr[0][alpha][gamma]+qr[1][alpha][gamma]+qr[2][alpha][gamma]
		if sumi<0.99 or sumi>1.01:
			print('alpha', alpha, 'gamma', gamma, ':', 'norm =', sumi)
