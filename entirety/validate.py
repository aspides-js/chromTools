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
from MACS3.IO.Parser import BEDParser, BEDPEParser

# ------------------------------------
# Misc function
# ------------------------------------

def args_validator( options ):
	''' 
	Validate command line arguments:
		chromhmm; files; outdir; genome
	'''

	## logging
	logging.basicConfig(level = 20, format='%(levelname)-5s @ %(asctime)s: %(message)s ',\
		datefmt='%a, %d %b %Y %H:%M:%S',\
		stream=sys.stderr, \
		filemode="w")

	options.error = logging.critical
	options.warn = logging.warning
	options.debug = logging.debug
	options.info = logging.info

	## chromhmm
	if os.path.exists( options.chromhmm + 'ChromHMM.jar'):
		options.info( 'Path to ChromHMM jarfile is: '+options.chromhmm )
		options.chromhmmJar=options.chromhmm + 'ChromHMM.jar'
	else:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), options.chromhmm +'ChromHMM.jar' )
	

	## genome
	try:
		options.genome = options.chromhmm+"CHROMSIZES/"+gnm[options.genome]
		options.info('Path to genome chromosome sizes is: %s' % options.genome)
	except:
		options.warn("Genome not recognised. Available shortcuts for genome chromosome sizes are: %s" % ", ".join(list(gnm.keys())))
		sys.exit(1)

	## files (path)

	for f in options.files:
		if not os.path.exists(f):
			raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f)


	## outdir
	options.metadir = options.outdir+"/0_metadata/"
	options.subdir = options.outdir+"/1_subsample/"
	options.bindir = options.outdir+"/2_binarised/"

	if options.force:
		try:
			os.makedirs(options.metadir, exist_ok = True)
			os.makedirs(options.subdir, exist_ok = True)
			os.makedirs(options.bindir, exist_ok = True)
		except:
			options.warn( "Output directories in %s could not be accessed. Terminating program." % options.outdir )
			sys.exit()
	else:
		try:
			os.makedirs(options.metadir)
			os.makedirs(options.subdir)
			os.makedirs(options.bindir)
		except:
			options.warn( "Output directories in %s could not be accessed or already exist. Please use --force-overwrite if you wish to overwrite output files. Terminating program." % options.outdir )
			sys.exit()		


	## seed
	if options.seed == None:
		options.seed=random.randint(0, 5000000000)
		options.info('RANDOM SEED: '+str(options.seed))
	else:
		options.info('RANDOM SEED: '+str(options.seed))

	## paired
	options.info("PE_MODE: %s" % options.paired)

	return options


def assert_compressed(f):
	try:
		open(f, 'r').readline()
		return False
	except:
		return True


def macs_validator( n, options ):
	''' 
	Set options for macs binarisation/peak calling
	'''

	## logging
	logging.getLogger().setLevel(30)

	options.error = logging.critical
	options.warn = logging.warning
	options.debug = logging.debug
	options.info = logging.info

	if (len(str(n)) < 2):
		nname = "0"+str(n)
	else:
		nname = str(n)


	## load_frag_files_options
	options.tfile = [options.subdir+'/downsampled.'+str(n)+'.bed']
	options.name = "P0"+nname
	options.cfile = False
	options.parser = BEDPEParser
	options.buffer_size = 100000

	## peakdetect options
	options.PE_MODE = options.paired
	options.log_qvalue = 5e-2
	options.log_pvalue = None
	options.maxgap = False
	options.minlen = False
	options.shift = False
	options.gsize = 2.9e9
	options.gsize = 10050000
	options.nolambda = False 
	options.smalllocal = 1000
	options.largelocal = 10000

	## call_peaks options
	options.store_bdg = False
	options.do_SPMR = False
	options.cutoff_analysis = False
	options.cutoff_analysis_file = "None"
	options.call_summits = False
	options.trackline = False
	options.broad = False

	# output filenames
	options.peakBed = os.path.join( options.outdir, options.name+"_peaks.bed" ) 
	options.bdg_treat = os.path.join( options.outdir, options.name+"_treat_pileup.bdg" ) 
	options.bdg_control= os.path.join( options.outdir, options.name+"_control_lambda.bdg" )

	options.fecutoff = 1.0

	return options

