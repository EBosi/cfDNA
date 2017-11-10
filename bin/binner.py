# Module 1: Make binning for a given window size

###########
import sys
from cPickle import load
import os
from IPython import embed
###########


###########
if __name__ == "__main__":
	usage = "python binner.py bin_size output [masked_repeats (T)]"
	# i.e. eb@eb-genetica-lab:~/projects/cfDNA$ python bin/binner.py 50000 bins_50K/bins.tsv
	args = sys.argv[1:]
	if len(args) == 2: bin_size,output,masking = args + [False]
	elif len(args) == 3: bin_size,output,masking = args + [True]
	else: exit("wrong usage ")
###########

## step 1) define chromosome numbers and partitions

def getChromoInfo(chromo_info="/home/eb/bioinfoSoftware/hg19/hg19_chromo_info.tsv"):
	""" parse a file containing chromosome name and length. Returns this info in a dict. """
	chromos = {}
	with open(chromo_info) as fh:
		for l in fh:
			k,v = l.strip().split()
			chromos[k] = int(v)
	return chromos

## step 2) produce a list of bins with the following formatting:
# chr; start; end; mapping;

def get_chromo_order():
	""" get an ordered iterations from a dictionary based on a filter """
	myOrder = map(lambda x: "chr%s" %x,range(1,23) + ["X","Y"])
	return myOrder

def get_masked_file(chromo,n_locationsDir="/home/eb/bioinfoSoftware/hg19/maskeds/"):
	""" get the file containing locations of Ns to define masked regions. """
	for f in os.listdir(n_locationsDir):
		if not f.endswith("n_locations"): continue 
		if f.split('.fa.')[0] == chromo: return os.path.abspath(os.path.join(n_locationsDir,f))

def binner(replicon,replicon_length,bin_size,masked_file="/home/eb/projects/cfDNA/gap.bed"):
	""" produce a list of start/end/masked tuples """
	from itertools import izip
	#if replicon == "chr17": embed()
	starts = range(1,replicon_length + 1,bin_size)
	ends = map(lambda x: x-1,starts[1:]) + [replicon_length]
	startEndGenerator = izip(starts,ends)
	start,end = startEndGenerator.next()
	yielded = set()
	with open(masked_file) as fh:
		for l in fh:
			l = l.strip().split()
			chr_ = l[1]
			if chr_ != replicon: continue
			s,e = map(lambda x: int(x) + 1,l[2:4])
			if e < start: continue
			while end < s:
				if start not in yielded: yield "%s\t%s\t%s\t0\n" %(replicon,start,end);yielded.add(start)
				try: start,end = startEndGenerator.next()
				except StopIteration: return
			if start < s < end:
				diff = end - s + 1
				if start not in yielded: yield "%s\t%s\t%s\t%s\n" %(replicon,start,end,diff);yielded.add(start)
				yielded.add(start)
				try: start,end = startEndGenerator.next()
				except StopIteration: return
				#continue			
			while e > start:
				if start < e <= end:
					diff = e - start
					if start not in yielded: yield "%s\t%s\t%s\t%s\n" %(replicon,start,end,diff);yielded.add(start)
					yielded.add(start)
					try: start,end = startEndGenerator.next()
					except StopIteration: return

				elif e >= end:
					diff = end - start + 1
					if start not in yielded: yield "%s\t%s\t%s\t%s\n" %(replicon,start,end,diff);yielded.add(start)
					yielded.add(start)
					try: start,end = startEndGenerator.next()
					except StopIteration: return
	while True:
		yield "%s\t%s\t%s\t0\n" %(replicon,start,end)
		try: start,end = startEndGenerator.next()
		except StopIteration: return
		

def getFinalTable(chromos,bin_size,out_table):
	""" write a table with the following formatting:
		chr; start; end; mapping;
		to the file `out_table` based on a given bin_size """
	with open(out_table,'w') as fh:
		for chromo in get_chromo_order():
			chr_size = chromos[chromo]
			print "working on %s ..." %chromo
			# masked_file = get_masked_file(chromo)
			for line in binner(chromo,chr_size,bin_size): fh.write(line)


########
# MAIN #
########

if __name__ == "__main__":
	bin_size = int(bin_size)
	chromo_d = getChromoInfo()
	getFinalTable(chromo_d,bin_size,output)


