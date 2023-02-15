#!/usr/bin/python

# ------------------------------------
# python modules
# ------------------------------------

import sys, os, time
import errno
import shutil
import random
import logging

# ------------------------------------
# own python modules
# ------------------------------------

from entirety.constants import GENOME as gnm

# ------------------------------------
# Misc function
# ------------------------------------

def args_validator( options ):
	''' 
	Validate command line arguments:
		chromhmm; files; outdir; genome
	'''

	## chromhmm
	if os.path.exists( options.chromhmm + 'ChromHMM.jar'):
		print( 'Path to ChromHMM jarfile is: '+options.chromhmm )
		options.chromhmmJar=options.chromhmm + 'ChromHMM.jar'
	else:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), options.chromhmm +'ChromHMM.jar' )
	
	## genome
	try:
		options.genome = options.chromhmm+"CHROMSIZES/"+gnm[options.genome]
		print('Path to genome chromosome sizes is: %s' % options.genome)
	except:
		print("Genome not recognised. Available shortcuts for genome chromosome sizes are: %s" % ", ".join(list(gnm.keys())))
		sys.exit(1)

	## files (path)

	for f in options.files:
		if not os.path.exists(f):
			raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f)


	## outdir
	options.metadir = options.outdir+"/0_metadata/"
	options.subdir = options.outdir+"/1_subsample/"
	options.bindir = options.outdir+"/2_binarised/"
	try:
		os.makedirs(options.metadir, exist_ok = True)
		os.makedirs(options.subdir, exist_ok = True)
		os.makedirs(options.bindir, exist_ok = True)
	except:
		sys.exit( "Output directories (%s) could not be accessed. Terminating program." % options.outdir )

	## seed
	if options.seed == None:
		options.seed=random.randint(0, 5000000000)
		print('RANDOM SEED: '+str(options.seed))
	else:
		print('RANDOM SEED: '+str(options.seed))

	return options

def assert_compressed(f):
	try:
		open(f, 'r').readline()
		return False
	except:
		return True