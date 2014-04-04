#! /usr/bin/env python
#copyright Simon Gravel, McGill University
#
#This is a bare bones version--if you would like to add features, let me know! 
#

#this python script loads in a site-frequency spectrum, and calculates the predicted number of nonreference or segregating sites to be discivered if a total of n samples are sequenced, using linear programming.

import sys
import numpy
from scipy.special import gammaln
from pulp import *
usage="lp.py file.txt -n totalsamples -b bins -B bootstraps -s solver [-o outputfile: default file.lp]  \n file.txt should be a site-frequency spectrum, a tab-delimited file showing the number of variants observed 0, 1, 2, ... , k times. Only the bins 2, 3,...,b+1 are used.\n totalsamples should be a list of comma separated ints that represent the target sizes. We assume that the first and last bins of the SFS are fixed sites. The actual values in these bins is not usedsolver is an external lp solver that python knows where to find. Currently only gurobi is implemented. Example usage ./lp.py samplesfs.txt -n 200,400 -B 100 -s gurobi -o samplesfs.lp"

args=sys.argv

if len(args)!=10 and len(args)!=12:
	print "found " ,len(args)-1, " arguments"
	print "incorrect number of fields"
	print usage
	sys.exit()

if args[2] !='-n' or args[4]!='-b' or args[6]!='-B' or args[8]!='-s' or (len(args)==12 and args[10]!='-o') :
	print "incorrect flags"
	
	print args
	print usage
	sys.exit()

seed=numpy.random.random_integers(1000000)
print "seed is", seed
numpy.random.seed(seed)

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
startcut=int(args[5])
nboot=int(args[7])
solver=args[9]
if(len(args)==10): 
	outf=args[9]
else:
	froot='.'.join(args[1].split('.')[:-1])
	outf=froot+".lp"

print "found a total of", numpy.sum(sfs), " sites,  including ", numpy.sum(sfs[1:-1]), " segregating, and ", sfs[0], " fixed ref, and " , sfs[-1], "fixed nonref" 
print "frequency spectrum contains ", len(sfs), "bins, assuming ", len(sfs)-1, " sequenced haploid genomes"	

n=len(sfs)-1

if min(Ns)<n:
	print "One of the extrapolation values is lower than the observed counts! You should be performing a projection rather than an extrapolation! (not yet implemented here)"
	sys.exit()

#projection matrix. Inspired by similar code in dadi by Gutenkunst (PLoS genetics, 2009). 

def lncomb(N,k):
    """
    Log of N choose k.
    """
    return gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1)

def projection(proj_to, proj_from, hits):
    """
    Coefficients for projection from a different fs size.

    proj_to: Number of samples to project down to.
    proj_from: Number of samples to project from.
    hits: Number of derived alleles projecting from.
    """
    if numpy.isscalar(proj_to) and numpy.isscalar(proj_from) and proj_from < proj_to:
        # Short-circuit calculation.
        print "error with projection values"
        raise
        contrib = numpy.zeros(proj_to+1)
    else:
        # We set numpy's error reporting so that it will ignore underflows, 
        # because those just imply that contrib is 0.
        previous_err_state = numpy.seterr(under='ignore', divide='raise',
                                          over='raise', invalid='raise')
        proj_hits = numpy.arange(proj_to+1) #array of target values
        # For large sample sizes, we need to do the calculation in logs, and it
        # is accurate enough for small sizes as well.
        lncontrib = lncomb(proj_to,proj_hits) #note: here we use the symmetry of the Hypergeometric distribution
        lncontrib += lncomb(proj_from-proj_to,hits-proj_hits)
        lncontrib -= lncomb(proj_from, hits)
        contrib = numpy.exp(lncontrib)
        numpy.seterr(**previous_err_state)
    return contrib
    
#define the projection matrix
def projmat(nto,nfrom):
	return numpy.array([projection(nto,nfrom,i) for i in range(nfrom+1)]).transpose()	


# Coarsen by keeping the first d bins. then use the following 2,4,8, etc 
def coarsenSFS(sfs, cut):
	#the SFS has a zeroton bin. We keep it for now, but will have to be removed later
	temp = list(sfs[:cut+1]); 
	
	maxi = cut + 1
	d = 1;
	curr = maxi; 
	while (maxi < len(sfs)):
		maxi += 2**d; 
		d +=1; 
		
   		temp.append(numpy.sum(sfs[(curr):(min(maxi, len(sfs)))]))
   		curr = maxi
   	return temp   

def coarsenSFS_straight(sfs, cut):
	#the SFS has a zeroton bin. We keep it for now, but will have to be removed later
	temp = list(sfs[:cut+1]); 
	temp.append(numpy.sum(sfs[cut+1:]))
   	return temp  

# Coarsen the projection matrixby keeping the first d bins. then use the following 2,4,8, etc
def coarsenMAT(mat, cut):
	#the matrix has no zeroton bin. 
	temp = list(mat[:cut,:]) 
	maxi = cut
	d = 1;
	curr = maxi; 
	while (maxi < len(mat)):
		maxi += 2**d; 
		d += 1; 
		
   		temp.append(numpy.sum(mat[(curr):(min(maxi, len(mat))),:],axis=0))
   		curr = maxi
   	return numpy.array(temp)
   	
def coarsenMAT_straight(mat, cut):
	#the matrix has no zeroton bin. 
	temp = list(mat[:cut,:]) 
	temp.append(numpy.sum(mat[cut:,:],axis=0))
   	return numpy.array(temp)

if solver=="gurobi":
	solvera=solvers.GUROBI_CMD(mip=False,msg=0)
	solverb=solvers.GUROBI_CMD(mip=False,msg=0)
else:
	print "error: solver ", solver, " not implemented"
	raise()
	
#a function that seeks the maximum number of bins that we can keep and still find a solution
def optimizecut(start,sfs,nmono,extrapto,bchp):
	keeptrying = True;
	step = 1; #step size by which to change the cutoff initially
 	extrapfrom=len(sfs)-1
	trycut = start + step;

	#bloc=projmat(extrapfrom,extrapto)#When extrapolating, the "to" and "from" are inverted
	#bchp=bloc[1:,1:]
	c=[1 for i in range(extrapto)]
	while keeptrying:
		trycut = trycut - step;
		if trycut == 0:
			break;
		subcoarse = coarsenSFS(sfs, trycut);
		obs = numpy.sum(subcoarse[1:]);
		bchpCoarse = coarsenMAT(bchp, trycut);
		nconst=bchpCoarse.shape[0]
		prob=LpProblem("maximum",LpMaximize)
		probmin=LpProblem("maximum",LpMinimize)
		#create the problem variables
		keys=map(lambda i:"%.3d" % i,range(1,extrapto+1))
		vmin=LpVariable.dict("phi",keys,0)
		vmax=LpVariable.dict("phi",keys,0)
   		#the SFS constraints
   		for i in range(nconst):
   			prob+=lpSum([bchpCoarse[i,j]*vmax[keys[j]] for j in range(extrapto)])==subcoarse[i+1]
   			probmin+=lpSum([bchpCoarse[i,j]*vmin[keys[j]] for j in range(extrapto)])==subcoarse[i+1]
		#the monotonicity constraints
   		for i in range(nmono):
   			prob+=vmax[keys[i+1]]>=vmax[keys[i+2]]
   			probmin+=vmin[keys[i+1]]>=vmin[keys[i+2]]
   		
   		#now define the target function
   		prob+=lpSum([c[j]*vmax[keys[j]] for j in range(extrapto)])
   		probmin+=lpSum([c[j]*vmin[keys[j]] for j in range(extrapto)])
   		#prob.solve
   		#print "Status:",LpStatus[prob.status] 
		
		ra=solvera.solve(prob)
		rb=solverb.solve(probmin)
   		if ra==1 and rb==1:#solution found!
   			keeptrying=False
   		
   		if not keeptrying:
			#print "cut",trycut,"range",value(probmin.objective),"-",value(prob.objective)
			return obs,value(probmin.objective),value(prob.objective),trycut

	


def extrapolate(sfs, popsize, startcut, nboot, propmon):
 
	sampsize = len(sfs) - 1;
	nmono = int(round(propmon*popsize));
	if (popsize == sampsize): #Then project down rather than extrapolate.
		print "Warning:Extrapolation size same as sample size"
   		tot = numpy.sum(sub[1:])
		return(tot);
  
  	if popsize < sampsize:  #Then project down rather than extrapolate.  
		print "Warning:Extrapolation size smaller than sample \
size:downsampling";
		bloc = projmat(sampsize, popsize);
		subsub = numpy.dot(bloc,sub); 
		tot = numpy.sum(subsub[1:]);
		return(tot);
  
  #Now that the pathological cases have ben dealt with, 
  #we can focus on the actual extrapolation problem*)
  
  #Projection matrices. Define here and pass faster than re-generate in the bootstrap
	bloc=projmat(sampsize,popsize)#When extrapolating, the "to" and "from" are inverted
	bchp=bloc[1:,1:]
  	ret1 = optimizecut(startcut,sfs,nmono,popsize,bchp)
	
	def fmap(i):
		poissonBoot = map(numpy.random.poisson, sfs);
		return optimizecut(startcut,poissonBoot,nmono,popsize,bchp)
	
	
	
	
	
	boots=map(fmap,range(nboot))
	return([ret1, boots])
  


fp=open(outf,"w")
fp.write("observed, min, max, number of bins used\n")
for N in Ns:
	point, boots=extrapolate(sfs, N, startcut, nboot, 0)
	print "estimate obtained"
	boots2=numpy.array([numpy.array(boot)*point[0]/boot[0] for boot in boots]) #normalize so that the initial number is the actual initial number observed. 
	#print boots, boots2
	
	CIs=numpy.percentile(boots2,[2.5,50,97.5],axis=0)


	print CIs
	
	fp.write("N=%d\n" % N)
	print point
	fp.write("\t".join(map(lambda i:"%f" % (i,), point))+"\n")
	fp.write("bootstrap, N=%d\n"% N)
	for boot in boots:
		fp.write("\t".join(map(lambda i:"%f" % (i,), boot))+"\n")
	fp.write("95% CI, \n")
	fp.write("\t".join(map(lambda i:"%f" % (i,), CIs[0]))+"\n")
	fp.write("\t".join(map(lambda i:"%f" % (i,), CIs[1]))+"\n")
	

fp.close()



