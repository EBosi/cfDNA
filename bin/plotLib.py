#########
import matplotlib.pyplot as plt
import numpy as np
#########

def boxPlot(d):
	""" wrap numpy's boxplot """
	xs,data = zip(*d.iteritems())
	plt.boxplot(data, positions=xs, notch=True)
	return


def scatterPlot(d):
	""" wrap plt scatterplot """
	xs = [i for i in d.keys() for j in d[i]]
	ys = [j for i in d.keys() for j in d[i]]
	plt.scatter(xs,ys)
	return 


def medianCorrection(d):
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
	d_out = {}
	m = np.median([j for i in d.values() for j in i])
	m_s = {k:np.median(v) for k,v in d.iteritems()}
	for k,v in d.iteritems():
		d_out[k] = map(lambda x: x*(m/float(m_s[k])),v)
	return d_out

def compareNormProcedures(roundiness=1,binSize=50000):
	import matplotlib.pyplot as plt
	import numpy as np
	import seaborn as sns
	f = "195-17.binned.50000.counts.reads.csv"
	binSize = float(binSize)
	gcs,counts,normalizeds = [],[],[]
	with open(f) as fh:
		for i,l in enumerate(fh):                        
			l = l.strip().split(',')                
			gc,count,norm = l[4],l[-4],l[-3]
			if l[5] == '0' and i!= 0: gcs.append(int(gc)); counts.append(float(count)); normalizeds.append(float(norm))

	d_gc,d_norm = {},{}
	xs,ys,norms = zip(*[(100*round(gc/binSize,roundiness),count,norm)
						for gc,count,norm in zip(gcs,counts,normalizeds)])
	for x,y in zip(xs,ys): d_gc[x] = d_gc.get(x,[]) + [y]
	d_norm = medianCorrection(d_gc)
	my_norms = [v for k in sorted(set(xs)) for v in d_norm[k]]
	
	data_y = list(ys) + list(norms) + list(my_norms)
	data_x = list(xs)*3
	
	hues = ["10M_norm","GC_lowess","GC_median"]
	hues = [h for h in hues for x in xs]
	
	sns.boxplot(x=data_x,y=data_y,hue=hues)


## ADD MEDIAN NORM


		
