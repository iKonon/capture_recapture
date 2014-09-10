#! /usr/bin/env python

#this python script loads in a site-frequency spectrum, and calculates the predicted number of segregating sites to be discivered if an additional n samples are sequenced, using a jackknife of order o.

import sys
import numpy


usage="jackknife.py file.txt -n totalsamples -p jackknifeorder [-o outputfile: default out.jk]  \n file.txt should be a site-frequency spectrum, a tab-delimited file showing the number of variants observed 0, 1, 2, ... , k times. Only the bins 2, 3,...,o+1 are used.\n totalsamples should be a list of comma separated ints. We assume that the first and last bins are not segregating. The actual values in these bins is not used, except as a santy check that the total number of sites predicted is not larger than the total number of sites available. example usage ./jackknife.py CEU.SFS -n 1000,2000,3000 -p 3 -o foo.jk"

args=sys.argv

if len(args)!=6 and len(args)!=8:
	print "found " ,len(args)-1, " arguments"
	print "incorrect number of fields"
	print usage
	sys.exit()
if args[2] !='-n' and args[4]!='-p' or (len(args)==8 and args[6]!='-o') :
	print "incorrect flags"
	print usage
	sys.exit()

try:
	fp=open(args[1],'r')
except IOError:
	print "could not open file ", args[1]
	sys.exit()
line=fp.readline()
lsp=line.split()
sfs= map(float,lsp)
Ns=args[3].split(',')
Ns=map(int,Ns)
order=int(args[5])
if(len(args)==8): 
	outf=args[7]
else:
	outf="out.jk"

print "found a total of", numpy.sum(sfs), " sites,  including ", numpy.sum(sfs[1:-1]), " segregating."
print "frequency spectrum contains ", len(sfs), "bins, assuming ", len(sfs)-1, " sequenced haploid genomes"	

n=len(sfs)-1

if min(Ns)<n:
	print "One of the extrapolation values is lower than the observed counts! You should be performing a projection rather than an extrapolation! (not yet implemented here)"
	sys.exit()

euler=0.5772156649
def harmonic(n):
	if(n<20):
		return numpy.sum(1./numpy.arange(1,n+1))
		
	else:
		return numpy.log(n)+euler+1./(2*n)-1./(12*n**2)+1./(120*n**4)

def deltaf(N,n):
	return harmonic(N-1)-harmonic(n-1)

outputs=[]
for N in Ns:
	delta=deltaf(N,n)
	
	if order==1:	
		missed=(-1.+n)/n*delta*sfs[1]
		
		
	if order==2:
		missed=(((1+2*(-2+n)*n)*delta)/(n*(-3+2*n))+((-2+n)*(-1+n)*delta**2)/(n*(-3+2*n)))*sfs[1]+(-((2*(-2+n)**2*delta)/((-1+n)*n*(-3+2*n)))-(2*(-2+n)**2*delta**2)/(n*(-3+2*n)))*sfs[2]
		
		
	if order==3:
		n=float(n)
		missed=(delta*(-127+628*n-1386*n**2+1722*n**3-1269*n**4+546*n**5-126*n**6+12*n**7+(-1+n)*(-10-87*n+310*n**2-372*n**3+210*n**4-57*n**5+6*n**6)*delta+2*(-2+n)**2*(-1+n)**3*(6-5*n+n**2)*delta**2)*sfs[1])/((-1+n)**2*n*(-5+2*n)*(-3+2*n)*(11-12*n+3*n**2))+(delta*(-2*(-248+402*n-61*n**2-249*n**3+195*n**4-57*n**5+6*n**6)-2*(-1+n)*(-560+1130*n-731*n**2+55*n**3+127*n**4-51*n**5+6*n**6)*delta+2*(-1+n)**2*(6-5*n+n**2)*(52-78*n+38*n**2-6*n**3)*delta**2)*sfs[2])/((-1+n)**2*n*(-5+2*n)*(-3+2*n)*(11-12*n+3*n**2))+(delta*(6*(3-2*n)**2*(-3+n)**3+6*(-3+n)**3*(-1+n)*(15-19*n+6*n**2)*delta+6*(-3+n)**2*(-1+n)**2*(-3+2*n)*(6-5*n+n**2)*delta**2)*sfs[3])/((-1+n)**2*n*(-5+2*n)*(-3+2*n)*(11-12*n+3*n**2))
	else:
		print "jackknife order ", o,"not implemented"
		sys.exit()
		
	print "based on %d-order jackknife, " % (order,), missed, "variants will be discovered if ", N,"samples are sequenced"	
	print "total nonreference will be", missed+numpy.sum(sfs[1:])
	outputs.append(missed+numpy.sum(sfs[1:]))
		
fp=open(outf,'w')		 
for i,j in zip(Ns,outputs):
	fp.write("%d\t%f\n"%(i,j))
fp.close()

