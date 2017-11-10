import pysam,sys
from IPython import embed

"""
samfile = pysam.AlignmentFile("1003-17.merged.bam", "rb")
myIter = samfile.fetch()
alls = {}
duplicated = {}

for i in myIter:
    if i.query_name in alls:
        old = alls[i.query_name]
        duplicated[i.query_name] = duplicated.get(i.query_name,[old]) + [i]
    alls[i.query_name] = i
   
for r in duplicated:
"""

def filteringCondition(alignment,MQ_threshold=15):
	""" filter alignment to get only unique reads """
	condition = alignment.flag / 256 == 0 # remove secondary alignment and other shit
	if alignment.mapping_quality < MQ_threshold: condition = False
	return condition
	
def getLeftMostPosition(alignment):
	""" returns 1-based leftmost position tuple (replicon,start) """
	return alignment.reference_start + 1

def createBinList(bin_file,indices):
	""" create a data structure comprising bins and counts (all 0 at the start)
		[[start1,end1,...,0],[start2,end2,...,0], # chr1
		 [start1,end1,...,0],[start2,end2,...,0], # chr2
		 ...]"""
	out = [[] for i in indices]
	with open(bin_file) as fh:
		for l in fh:
			l_ = l.strip().split()
			out[indices[l_[0]]].append(l_ + [0])
	return out

def getBin(start,binSize,binList):
	""" get bin corresponding to a given position """
	indx,left = start/binSize , start%binSize
	try: binList[indx]
	except IndexError:
		if indx == len(binList): return None
		else: raise IndexError
	if left != 0: return binList[indx]
	return binList[indx+1]

def addToBin(alignment,binSize,binList,binIndx):
	""" check for conditions and add to binList """
	if not filteringCondition(alignment): return
	position = getLeftMostPosition(alignment)
	sublist = binList[binIndx]
	try: bin_ = getBin(position,binSize,sublist)
	except: embed()
	if bin_: bin_[-1] += 1

def main(bin_file,bam,bin_size):
	my_bam = pysam.AlignmentFile(bam, "rb")
	chromo2Index = {v:i for i,v in enumerate(map(lambda x: 'chr%s' %x,range(1,23)) + ['chrX','chrY'])}
	#embed()
	#sys.exit()
	bin_list = createBinList(bin_file,chromo2Index)
	for al in my_bam.fetch():
		binIndx = chromo2Index.get(al.reference_name)
		if binIndx == None: continue
		addToBin(al,bin_size,bin_list,binIndx)
	return bin_list

if __name__ == "__main__":
	usage = "python bam_analyser.py binFile bamFile binSize outFile"
	try: bin_file,bam,bin_size,outfile = sys.argv[1:]
	except: print usage; sys.exit()
	bin_size = int(bin_size)
	bin_list = main(bin_file,bam,bin_size)
	with open(outfile,'w') as fh:
		for l in bin_list:
			for l2 in l:
				l2[-1] = str(l2[-1])
				fh.write('\t'.join(l2)+'\n')
	
	
