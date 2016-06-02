#! /usr/bin/env python

#this python script loads in a site-frequency spectrum, and calculates the predicted number of segregating sites to be discivered if an additional n samples are sequenced, using a jackknife of order o.

import sys
import numpy


usage = "jack.py file.txt -n totalsamples -p jackknifeorder [-o outputfile: default out.jk]\
  \n file.txt should be a site-frequency spectrum, a tab-delimited file showing the number\
   of variants observed 0, 1, 2, ... , k times. Only the bins 1, 2,...,p are used.\n\
    totalsamples should be a list of comma separated ints. We assume that the first and\
    last bins are not segregating. The actual values in these bins is not used, except as\
     a santy check that the total number of sites predicted is not larger than the total\
      number of sites available. example usage ./jackknife.py CEU.SFS -n 1000,2000,3000\
       -p 3 -o foo.jk"

args = sys.argv
print(args)
if len(args) != 6 and len(args) != 8:
	print "found " , len(args) - 1, " arguments"
	print "incorrect number of fields"
	print usage
	sys.exit()
if args[2] != '-n' and args[4] != '-p' or (len(args) == 8 and args[6] != '-o'):
	print "incorrect flags"
	print usage
	sys.exit()

try:
	fp = open(args[1], 'r')
except IOError:
	print "could not open file ", args[1]
	sys.exit()
line = fp.readline()
lsp = line.split()
sfs = map(float,lsp)
Ns = args[3].split(',')
Ns = map(int, Ns)
order = int(args[5])
if(len(args) == 8): 
	outf = args[7]
else:
	outf = "out.jk"

print "found a total of", numpy.sum(sfs), " sites,  including ", numpy.sum(sfs[1:-1]), " segregating."
print "frequency spectrum contains ", len(sfs), "bins, assuming ", len(sfs)-1, " sequenced haploid genomes"	

n = len(sfs)-1

if min(Ns) < n:
	print "One of the extrapolation values is lower than the observed counts! You should be performing a projection rather than an extrapolation! (not yet implemented here)"
	sys.exit()

euler = 0.5772156649
def harmonic(n):
	if(n < 20):
		return numpy.sum(1. / numpy.arange(1, n+1))
		
	else:
		return numpy.log(n) + euler + 1./(2*n) - 1./(12*n**2) + 1./(120*n**4)

def deltaf(N,n):
	return harmonic(N-1) - harmonic(n-1)

outputs = []
for N in Ns:
	delta = deltaf(N,n)
	
	if order == 1:	
		missed = (-1.+n) / n * delta * sfs[1]
		
		
	elif order == 2:
		missed = (((1+2*(-2+n)*n)*delta)/(n*(-3+2*n))+((-2+n)*(-1+n)*delta**2)/(n*(-3+2*n)))*sfs[1]+(-((2*(-2+n)**2*delta)/((-1+n)*n*(-3+2*n)))-(2*(-2+n)**2*delta**2)/(n*(-3+2*n)))*sfs[2]
		
		
	elif order == 3:
		n = float(n)
		missed = (delta*(-127+628*n-1386*n**2+1722*n**3-1269*n**4+546*n**5-126*n**6+12*n**7+(-1+n)*(-10-87*n+310*n**2-372*n**3+210*n**4-57*n**5+6*n**6)*delta+2*(-2+n)**2*(-1+n)**3*(6-5*n+n**2)*delta**2)*sfs[1])/((-1+n)**2*n*(-5+2*n)*(-3+2*n)*(11-12*n+3*n**2))+(delta*(-2*(-248+402*n-61*n**2-249*n**3+195*n**4-57*n**5+6*n**6)-2*(-1+n)*(-560+1130*n-731*n**2+55*n**3+127*n**4-51*n**5+6*n**6)*delta+2*(-1+n)**2*(6-5*n+n**2)*(52-78*n+38*n**2-6*n**3)*delta**2)*sfs[2])/((-1+n)**2*n*(-5+2*n)*(-3+2*n)*(11-12*n+3*n**2))+(delta*(6*(3-2*n)**2*(-3+n)**3+6*(-3+n)**3*(-1+n)*(15-19*n+6*n**2)*delta+6*(-3+n)**2*(-1+n)**2*(-3+2*n)*(6-5*n+n**2)*delta**2)*sfs[3])/((-1+n)**2*n*(-5+2*n)*(-3+2*n)*(11-12*n+3*n**2))
	elif order == 4:
		n = float(n)
		# The following expression was hand-converted from mathematica notebook 
		# analytical_checks_fourthorder.
		missed = (delta*(2*(79138 - 768163*n + 3305209*n**2 - 8455556*n**3 + 14456396*n**4\
		 - 17528889*n**5 + 15570611*n**6 - 10302336*n**7 + 5104884*n**8 - 1885839*n**9 \
		 + 511449*n**10 - 98781*n**11 + 12849*n**12 - 1008*n**13 + 36*n**14) \
		 + (-1 + n)*(-142984 + 984880*n - 3219708*n**2 + 6550921*n**3 - 9150649*n**4 \
		 + 9174455*n**5 - 6742096*n**6 + 3659027*n**7 - 1462263*n**8 + 424335*n**9 \
		 - 86856*n**10 + 11874*n**11 - 972*n**12 + 36*n**13)*delta + 2*(-1 + n)**2\
		 *(6 - 5*n + n**2)*(-2288 + 3592*n + 10090*n**2 - 36645*n**3 + 50160*n**4 \
		 - 39031*n**5 + 19017*n**6 - 5922*n**7 + 1147*n**8 - 126*n**9 + 6*n**10)*delta**2\
		 + (-4 + n)*(-1 + n)**3*(6 - 5*n + n**2)**2*(184 - 606*n + 803*n**2 - 549*n**3 \
		 + 204*n**4 - 39*n**5 + 3*n**6)*delta**3)*sfs[1])/((-2 + n)*(-1 + n)**3*n\
		 *(-7 + 2*n)*(-5 + 2*n)*(-3 + 2*n)*(5 - 5*n + n**2)*(26 - 18*n + 3*n**2)*(11\
		 - 12*n + 3*n**2)) + (delta*(2*(806336 - 4114800*n + 9788448*n**2 - 14623884*n**3\
		 + 15584978*n**4 - 12651005*n**5 + 8033192*n**6 - 3985285*n**7 + 1516278*n**8\
		 - 429717*n**9 + 87120*n**10 - 11877*n**11 + 972*n**12 - 36*n**13)\
		 - 2*(-1 + n)*(-1762976 + 8275744*n - 17540568*n**2 + 22472256*n**3\
		 - 19868590*n**4 + 13296025*n**5 - 7229740*n**6 + 3298055*n**7 - 1236311*n**8\
		 + 360392*n**9 - 76491*n**10 + 10944*n**11 - 936*n**12 + 36*n**13)*delta\
	     - 4*(-1 + n)**2*(6 - 5*n + n**2)*(-117832 + 426816*n - 645654*n**2 + 513018*n**3\
	     - 207385*n**4 + 15319*n**5 + 25646*n**6 - 13195*n**7 + 3075*n**8 - 366*n**9\
	    + 18*n**10)*delta**2 - 2*(-4 + n)*(-1 + n)**3*(6 - 5*n + n**2)**2*(3176 - 8768*n\
	    + 9854*n**2 - 5764*n**3 + 1850*n**4 - 309*n**5 + 21*n**6)*delta**3)*sfs[2])\
	    /((-2 + n)*(-1 + n)**3*n*(-7 + 2*n)*(-5 + 2*n)*(-3 + 2*n)*(5 - 5*n + n**2)\
	    *(26 - 18*n + 3*n**2)*(11 - 12*n + 3*n**2)) + (delta*(12*(3 - 2*n)**2*(-3 + n)\
	    *(24678 - 73356*n + 88865*n**2 - 54881*n**3 + 16400*n**4 - 464*n**5 - 1289*n**6\
	    + 419*n**7 - 57*n**8 + 3*n**9) + 6*(-1 + n)*(9 - 9*n + 2*n**2)\
	    *(-354348 + 1276800*n - 1961648*n**2 + 1662242*n**3 - 829199*n**4 + 230491*n**5\
	     - 21616*n**6 - 6925*n**7 + 2631*n**8 - 354*n**9 + 18*n**10)*delta +6*(-1 + n)**2\
	     *(-3 + 2*n)*(6 - 5*n + n**2)*(138996 - 425226*n + 547888*n**2 - 385039*n**3\
	     + 158306*n**4 - 36911*n**5 + 3767*n**6 + 208*n**7 - 87*n**8 + 6*n**9)*delta**2
	     + 6*(-4 + n)*(-1 + n)**3*(-3 + 2*n)*(6 - 5*n + n**2)**2*(-1494 + 2874*n\
	     - 2121*n**2 + 758*n**3 - 132*n**4 + 9*n**5)*delta**3)*sfs[3])/((-2 + n)\
	     *(-1 + n)**3*n*(-7 + 2*n)*(-5 + 2*n)*(-3 + 2*n)*(5 - 5*n + n**2)*(26 - 18*n\
	     + 3*n**2)*(11 - 12*n + 3*n**2)) + (delta*(-12*(3 - 2*n)**2*(-4 + n)**4*(-3 + n)\
	     *(11 - 12*n + 3*n**2)**2 - 12*(-4 + n)**4*(-1 + n)*(9 - 9*n + 2*n**2)\
	     *(-803 + 2196*n - 2363*n**2 + 1249*n**3 - 324*n**4 + 33*n**5)*delta\
	     - 24*(-4 + n)**4*(-1 + n)**2*(-3 + 2*n)*(6 - 5*n + n**2)*(143 - 299*n + 228*n**2\
	     - 75*n**3 + 9*n**4)*delta**2 - 12*(-4 + n)**4*(-1 + n)**3*(-3 + 2*n)\
	     *(6 - 5*n + n**2)**2*(11 - 12*n + 3*n**2)*delta**3)*sfs[4])/((-2 + n)\
	     *(-1 + n)**3*n*(-7 + 2*n)*(-5 + 2*n)*(-3 + 2*n)*(5 - 5*n + n**2)\
	     *(26 - 18*n + 3*n**2)*(11 - 12*n + 3*n**2))
	else:
		print "jackknife order ", order, "not implemented"
		sys.exit()
		
	print "based on order %d jackknife, " % (order,), missed, " additional variants will be discovered if ", N, "samples are sequenced"	
	print "total nonreference will be", missed + numpy.sum(sfs[1:])
	outputs.append(missed+numpy.sum(sfs[1:]))
		
fp = open(outf,'w')		 
for i, j in zip(Ns,outputs):
	fp.write("%d\t%f\n" % (i, j))
fp.close()

