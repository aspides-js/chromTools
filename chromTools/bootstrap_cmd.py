#!/usr/bin/python


# ------------------------------------
# python modules
# ------------------------------------

import sys, os, time
import shutil
import mmh3
import logging
import pandas as pd
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt

# ------------------------------------
# Main function
# ------------------------------------
def run_plot( options):
	options.info('calualte n plot')
	plot_metrics()

def run_split(options):
	start_time = time.time()
	options.start_time = start_time

	options.nfile = 4


	pool = mp.Pool()
	args = [(n, options) for n in range(1, options.nfile)] 
	#pool.starmap(split, args)
	pool.starmap(learn_model, args)
	pool.close()
	options.info("--- %s seconds ---" % (time.time() - start_time))

	print('Complete')

	
def run_bootstrap_chmm(options):
	start_time = time.time()
	options.start_time = start_time

	options.nfile = 4

	options.bootdir = "/lustre/home/jms260/package/tmp/chip/split/"
	options.datatype = "h3k27ac"

	if not os.path.exists(options.subdir+"/downsampled.0.bed"):
		options.warning("Concatenated file not found so concatenating those in --files to %s" % options.bootdir)
		cat_bed( options.files, options.bootdir, options.info)

	pool = mp.Pool()
	args = [(n, options) for n in range(1, options.nfile)] 
	#pool.starmap(split, args)
	pool.starmap(learn_model, args)
	pool.close()
	options.info("--- %s seconds ---" % (time.time() - start_time))

#--------------------------------------------------------------------------------#

def split(n, options):
	''' 
	Split a bed file. For each read pair, a random value is assigned between the range. 
	The proportion is used to calculate a maximum acceptable value within the range. Records 
	whose value is below the limit are written to split1, records whose hash value is above 
	the limit are written to split2.

	'''
	options = assign_dirs(n, options)
	os.makedirs(options.n1dir, exist_ok = True)
	os.makedirs(options.n2dir, exist_ok = True)

	split1 = options.n1dir+'/split.bed'
	split2 = options.n2dir+'/split.bed'
	a=params(0.5)
	with open(split1, 'w') as out1f:
		with open(split2, 'w') as out2f:
			with open(options.subdir+"downsampled.0.bed",'r') as f:
				for line in f:
					if (discard(a, options.seed, line)):
						out1f.write(line)
					else:
						out2f.write(line)


def makeCHMMfile(options):
	os.makedirs(options.n1metadir, exist_ok = True)
	os.makedirs(options.n2metadir, exist_ok = True)

	dirDict = {options.n1metadir: options.n1dir, options.n2metadir: options.n2dir}

	for metadir, dir in dirDict.items():
		with open( metadir+'cellMarkFile.txt', 'w') as f:
			f.write("N+\t"+options.datatype+"\tsplit.bed")


def binarise(n, options):
	'''
	Binarise the data
	'''
	#split(n, options)
	options = assign_dirs(n, options)
	os.makedirs(options.n1bindir, exist_ok = True)
	os.makedirs(options.n2bindir, exist_ok = True)

	makeCHMMfile(options)
	
	dirDict = {options.n1metadir: [options.n1dir, options.n1bindir], options.n2metadir: [options.n2dir, options.n2bindir]}

	for metadir, dir in dirDict.items():
		print("java -mx2400M -jar %s BinarizeBed -b 200 %s %s %scellMarkFile.txt %s" % (options.chromhmmJar, options.genome, dir[0], metadir, dir[1]+options.datatype))
		#os.system("java -mx2400M -jar %s BinarizeBed -b 200 %s %s %scellMarkFile.txt %s" % (options.chromhmmJar, options.genome, dir[0], metadir, dir[1]+options.datatype))

	merge_binary(n, options)

def merge_binary(n, options):
	options.info("Merging different data modalities to single chromosome binary files")
	options = assign_dirs(n, options)

	dirDict = {options.n1metadir: [options.n1dir, options.n1bindir, options.n1mergedir], 
				options.n2metadir: [options.n2dir, options.n2bindir, options.n2mergedir]}

	for dir in dirDict.values():
		os.system("java -mx2400M -jar %s MergeBinary %s %s" % (options.chromhmmJar, dir[1], dir[2]))


def learn_model(n, options):
	options.info("Learning models")
	options = assign_dirs(n, options)

	dirDict = {options.n1metadir: [options.n1dir, options.n1bindir, options.n1mergedir, options.n1moddir], 
				options.n2metadir: [options.n2dir, options.n2bindir, options.n2mergedir, options.n2moddir]}

	for dir in dirDict.values():
		for state in range(options.state1, options.staten1):
			os.system("java -mx2400M -jar %s LearnModel -p 0 %s %s %s hg38" % (options.chromhmmJar, dir[2], dir[3], state))




#--------------------------------------------------------------------------------#

def read_emission_table(tablePath):
    t = pd.read_csv(tablePath, sep='\t')
    t=t.reindex(sorted(t.columns[1:]), axis = 1)
    t = (t>=0.5).astype(int)
    t = np.array(t)
    return t

def calc_corr(boot):
    boots = [boot,boot+1]
    metric=[]
    for state in range(5, 16):
        paths = ['/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/bootstrap/'+str(i)+'/3_model/emissions_'+str(state)+'.txt' for i in boots]
        emis = [read_emission_table(n) for n in paths]
        compare = [[np.any(np.all(x == emis[0], axis=1)) for x in emis[1]], 
                [np.any(np.all(x == emis[1], axis=1)) for x in emis[0]]]
        metric.append((1-((len(compare[0])*2)-sum(compare[0])-sum(compare[1]))/10))
    return metric 

def calc_euc(boot):
    boots = [boot,boot+1]
    euclidean = {}
    for state in range(5, 16):
        paths = ['/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/bootstrap/'+str(i)+'/3_model/emissions_'+str(state)+'.txt' for i in boots]
        emis = [read_emission_table(n) for n in paths]
        euc_all = []
        for x in range(0,len(emis[0])):
            euc = []
            for y in range(0,len(emis[1])):
                euc.append(np.linalg.norm(emis[0][x]- emis[1][y]))
            euc = min(euc/max(euc))
            euc_all.append(euc)
        euclidean.setdefault(len(emis[0]), [])
        euclidean[len(emis[0])].append(sum(euc_all)/len(euc_all))
    return euclidean

def plot_metrics():
	m=[]
	e = pd.DataFrame()
	for i in range(0, 20, 2):
		m.append(calc_corr(i))
		euc = pd.DataFrame.from_dict(calc_euc(i))
		e = pd.concat([e, euc])
		c = pd.DataFrame(m)
		c.columns = [*range(5,16)]
		c.boxplot()
		plt.savefig('/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/bootstrap/corrplot.jpg')
		plt.close()
		e.boxplot()
		plt.savefig('/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/bootstrap/eucplot.jpg')
