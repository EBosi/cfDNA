# normalizer

###########
import sys
from cPickle import load
import os
import numpy as np
import scipy as sp
from IPython import embed
###########

def getNormalizedCounts(f,autosome=True):
	""" get 10M normalized counts """
	# initialize variables
	inputs,finals = [[],[]],[]
	auto_inp,auto_final = [],[]
	XY_inp,XY_final = [],[]
	# parse file, populate inputs
	with open(f) as fh:
		for l in fh:
			l_ = l.strip().split()
			if l_[0] in {"chrX","chrY"}: inputs[1].append(l_)
			else: inputs[0].append(l_)
	# get normalized counts, compute GC
	if autosome: finals = map(normalize10M,inputs)
	else: finals = [normalize10M(inputs[0] + inputs[1])]
	# FIN
	return finals

def GC_analysis(normalized):
	""" wrap GC analysis for a list of normalized stuff """
	GC_distr = {}
	for n in normalized:
		if n[3] != '0': continue
		GC_frac = n[-1]
		GC_distr[GC_frac] = GC_distr.get(GC_frac,[]) + [n[-2]]
	return GC_distr

def countFile2GCDistr(count_list,autosome=True):
	""" parse file and produce a dict of GC distributions """
	normalized = normalize10M(count_list)
	GC_distribution = map(GC_analysis,normalized)
	return GC_distribution
	
def lowess(x, y, f=2. / 3., iter=3):
	from math import ceil
	import numpy as np
	from scipy import linalg

	"""lowess(x, y, f=2./3., iter=3) -> yest
	Lowess smoother: Robust locally weighted regression.
	The lowess function fits a nonparametric regression curve to a scatterplot.
	The arrays x and y contain an equal number of elements; each pair
	(x[i], y[i]) defines a data point in the scatterplot. The function returns
	the estimated (smooth) values of y.
	The smoothing span is given by f. A larger value for f will result in a
	smoother curve. The number of robustifying iterations is given by iter. The
	function will run faster with a smaller number of iterations.
	"""
	n = len(x)
	r = int(ceil(f * n))
	h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
	w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
	w = (1 - w ** 3) ** 3
	yest = np.zeros(n)
	delta = np.ones(n)
	for iteration in range(iter):
		for i in range(n):
			weights = delta * w[:, i]
			b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
			A = np.array([[np.sum(weights), np.sum(weights * x)],
						  [np.sum(weights * x), np.sum(weights * x * x)]])
			beta = linalg.solve(A, b)
			yest[i] = beta[0] + beta[1] * x[i]

		residuals = y - yest
		s = np.median(np.abs(residuals))
		delta = np.clip(residuals / (6.0 * s), -1, 1)
		delta = (1 - delta ** 2) ** 2
	return yest

def lowessCorrection(count,lowess_mean,exp_count):
	""" corrected number of reads using lowess approach.
		It is estimated as follows:
		r_c = r_u - (r_loess - e_r)
		where:
			r_c = corrected number of reads
			r_u = uncorrected number of reads
			r_loess = mean of read count corresponding to that GC bin smoothed with lowess
			e_r = expected count (average count for that chromosome) """
	# TODO
	return
	
def medianCorrection(l,d):
	""" compute a correction coefficient (cc) based on medians, m/m_x where:

		m = median of counts for all bins
		m_x = median of counts for a given bin (x)

		corrected read counts is count * cc
		
		requires:
			- l: list of count rows
			- d: dict of GC bins
		returns:
			- list of count rows with the corrected values as additional column
	 """
	m = np.median([j for i in d.values() for j in i])
	m_s = {k:np.median(v) for k,v in d.iteritems()}
	for i,l_ in enumerate(l):
		count,gc = l_[-2:]
		cc = m/float(m_s[gc])
		corrected_count = count * cc
		l_ += [corrected_count]
	return l
	
def GC_norm_lowess(d):
	""" wrapper for the `lowess` fxn """
	xs,ys = [],[]
	for k,v in d.iteritems():
		xs += [k] * len(v)
		ys += v
	yest = lowess(xs,ys)

###################################

def GC_fraction(row):
	""" compute the GC fraction for a bin """
	if isinstance(row,str): row = row.split()
	start,end,GC = map(float,row[1:3] + [row[-3]])
	return round(GC / (end - start + 1),1)*100

def normalize10M(l):
	""" returns a 10M normalized vector with GC """
	norm_factor = 10000000./sum([float(i[-1]) for i in l])
	for i in l:
		i += [int(i[-1]) * norm_factor]
		i[-1] *= norm_factor
		i += [GC_fraction(i)]
	return l

def setDataType(row):
	""" set all members but the first element as float """
	return [row[0]] + map(float,row[1:])

def GC_analysis(normalized):
	""" wrap GC analysis for a list of normalized stuff """
	GC_distr = {}
	for n in normalized:
		if n[3] != 0: continue
		GC_frac = n[-1]
		GC_distr[GC_frac] = GC_distr.get(GC_frac,[]) + [n[-2]]
	return GC_distr

def medianCorrection(l,d):
	""" compute a correction coefficient (cc) based on medians, m/m_x where:

		m = median of counts for all bins
		m_x = median of counts for a given bin (x)

		corrected read counts is count * cc
		
		requires:
			- l: list of count rows
			- d: dict of GC bins
		returns:
			- list of count rows with the corrected values as additional column
	 """
	m = np.median([j for i in d.values() for j in i])
	m_s = {k:np.median(v) for k,v in d.iteritems()}
	for i,l_ in enumerate(l):
		if l_[3] != 0: l_ += [0]; continue
		count,gc = l_[-2:]
		try: cc = m/float(m_s[gc])
		except: embed()
		corrected_count = count * cc
		l_ += [corrected_count]
	return l

def normalizeSingleFile(f,fxn=medianCorrection):
	""" wrapper to normalize single file """
	lines = [setDataType(i.strip().split()) for i in open(f)]
	normalized_10M = normalize10M(lines)
	gc_bins = GC_analysis(normalized_10M)
	return medianCorrection(normalized_10M,gc_bins)

def main(f,out):
	normalized = normalizeSingleFile(f)
	o = open(out,'w')
	o.write("Chromosome\tStart\tEnd\tNs\tGC\tRawCounts\tCounts10M\tCounts10MGC\n")
	for l in normalized:
		l_ = map(str,l)
		l_ = l_[:-2] + [l_[-1]]
		o.write("\t".join(l_) + "\n")

if __name__ == "__main__":
	inp,out = sys.argv[1:]
	main(inp,out)











