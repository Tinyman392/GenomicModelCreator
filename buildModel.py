#!/usr/bin/env python

'''
python buildModel.py [-<-> a<rg> val]
'''

from sys import stderr
import sys
from optparse import OptionParser
import os
import KMC
import Tabular
import PATRICPublicTabular
import xgboost as xgb
import shutil
import Matrix
import LibSVM
import TrainXGBoost
from ast import literal_eval
import ModelCleanup
import ModelMetrics
import TrainSciKit
import PLF
import Alignment
import glob

os.environ["CUDA_VISIBLE_DEVICES"]=""

def strToBool(s):
	if s.lower()[0] == 't':
		s = True
	else:
		s = False

	return s

def makeDir(d, rm = True):
	flag = True
	if not os.path.isdir(d):
		if os.path.exists(d):
			print d, 'is not a directory'
			flag = False
		else:
			os.mkdir(d)
	else:
		if rm:
			shutil.rmtree(d)
			os.mkdir(d)

	return flag

def checkDir(d):
	if os.path.isdir(d):
		return True
	else:
		print 'directory does not exist:', d
		return False

def checkFile(f):
	if os.path.exists(f):
		return True
	else:
		print 'file does not exist:', f
		return False

def cleanDir(d):
	if d == '':
		return d

	if d[-1] != '/':
		d += '/'

	return d

def getOptions():
	parser = OptionParser()

	parser.add_option('-f', '--fasta_dir', help="Directory name containing fasta files to train on.  File names must be [GID].fasta where [GID] matches genome IDs found in the first column of the tabular file or public genome database", metavar='DIR', default='', dest='fastaDir')
	parser.add_option('-t', '--tabular_file', help="Tabular file to use.  For private genome models, it's 3 columns: genome id, antibiotic:test_method, label.  For public genomes, use the PATRIC AMR Tabular file format", metavar='FILE', default='', dest='tabFile')
	parser.add_option('-H', '--header', help = 'Specify whether there is a header or not, defaults to False, there is no header', metavar='BOOL=False', default="False", dest='header')
	parser.add_option('-T', '--temp_dir', help='Directory name to store temporary files.  Defaults to "temp"', metavar='DIR=temp', default='temp', dest='tempDir')
	parser.add_option('-o', '--out_dir', help='Directory name to output files to.  Defaults to "model"', metavar='DIR=model', default='model', dest='outDir')
	parser.add_option('-n', '--threads', help="Number of threads to run.  Defaults to 1", metavar='INT=1', type=int, default=1, dest='threads')
	parser.add_option('-d', '--depth', help="Specify the depth of the XGBoost model, defaults to 4", metavar='INT=4', type=int, default=4, dest='depth')
	parser.add_option('-k', '--kmer_size', help="Specify the kmer size to use.  Defaults to 10", metavar='INT=10', type=int, default=10, dest='kmerSize')
	parser.add_option('-K', '--kmc_dir', help='Specify the outupt directory for pre-computed KMC output.  If supplied, KMC will not be run and this directory used instead.', metavar='DIR', default='', dest='kmcDir')
	parser.add_option('-p', '--public', help='Boolean value (true or false) to specify whether or not to use public genome data.  Defaults to "false"', metavar='BOOL=false', default='false', dest='public')
	parser.add_option('-g', '--genus', help='Specify the genus to train on for a public data model.  Separate multiple genera with commas.', metavar='GENUS', default='', dest='genus')
	parser.add_option('-s', '--test_standard', help='Specify the breakpoints testing standard to accept without conversion from the public genome data.  If multiple, separate with commas.  If none supplied, all taken.', metavar='TEST_STANDARD', default='', dest='testStandard')
	parser.add_option('-b', '--breakpoints', help="Specify the file containing the breakpoint calls.  It is 4 columns: taxonomy, antibiotic, S, R, I.", metavar='FILE', default='', dest='breakpointsFile')
	parser.add_option('-P', '--presence_absence', help="Specify whether or not to use presence vs absence of the kmer vs kmer counts", metavar='BOOL=false', default='false', dest='presence_absence')
	parser.add_option('-i', '--individual', help="Specify whether to run individual models per antibiotic", metavar='BOOL=false', default='false', dest='individualModels')
	parser.add_option('-F', '--filter_class', help="Specify minimum number of genomes required in 2 bins to run antibiotic model", type=int, metavar="INT=0", default=0, dest='filterClass')
	parser.add_option('-e', '--enumerate_classes', help="Specify whether or not to enumerate classes.  DO NOT ENUMERATE SIR LABELS!!!", default='False', dest='enumerateClasses')
	parser.add_option('-a', '--folds_to_run', help="Specify the number of folds to run.  Must be <= total folds (-A|--total_folds; default = 10).  Defaults to 5.", metavar='INT=5', type=int, default=5, dest='foldsToRun')
	parser.add_option('-A', '--total_folds', help='Specify the total folds to create.  Must be >= folds to run (-a|--folds_to_run; default = 5).  Defaults to 10.', metavar='INT=10', type = int, default=10, dest='totalFolds')
	parser.add_option('-c', '--classification', help="Specify whether or not to run a classification model.  Defaults to false", metavar='BOOL=false', default='false', dest='classify')
	parser.add_option('-j', '--SvNS', help = 'Specify to train a model of susceptible vs non-susceptible', metavar='BOOL=false', default='false', dest = 'svns')
	parser.add_option('-J', '--noI', help = 'Specify to train a model where the intermediate class is ignored', metavar='BOOL=false', default='false', dest = 'noI')
	parser.add_option('-m', '--model_params', help="Specify the parameters for the XGBoost model", metavar='PYTHON_HASH', default="{'eta':0.0625, 'silent':1}", dest='parameters')
	parser.add_option('-N', '--num_rounds', help="Specify the number of rounds to boost", metavar='INT=1000', type=int, default=1000, dest='num_rounds')
	parser.add_option('-C', '--cleanup', help='Specify whether to remove large files and directories after running', metavar='BOOL=true', default='true', dest='cleanup')
	parser.add_option('-S', '--compute_stats', help='Specify which statistics to compute on the model and output.  We currently support rawAcc, w1Acc, r2, VMEME, confMatrix, and classReport.  Combinations for AMRcls, AMRreg, and cls are supported as well.  Future support for heatmap will be added at a later date.  Separate multiple with commas.', metavar='LIST=""', default='', dest='stats')
	# parser.add_option('-u', '--shuffle', help='Specify whether to shuffle by genome ID (gid), nothing (all), or not shuffle (false).', metavar='gid|all|none=none', default='none', dest='shuffle')
	parser.add_option('-M', '--model_type', help="Specify the model type to train with: XGBoost, RandomForest, ExtraTrees, Bagging, or SVM.", metavar="XGBoost|RandomForest|ExtraTrees|Bagging|SVM=XGBoost", default="XGBoost", dest="modelType")
	parser.add_option('-E', '--max_features', help='Specify the maximum number of features for non-XGBoost ensemble methods', metavar='NUM=0.75', default=0.75, type=float, dest='max_features')
	parser.add_option('-l', '--max_raw_sample', help='Specify the number of rows to sample for non-XGBoost ensemble methods', metavar='NUM=0.75', default=0.75, type=float, dest='row_sample')
	parser.add_option('-O', '--stats_only', help='Specify whether to only run statistics on an already trained model', metavar='BOOL=false', default='false', dest='statsOnly')
	parser.add_option('-q', '--paired_end', help='Specify whether to use paired end reads as input.  Fasta file input should be a file containing a list of paired end read files for each genome.', metavar='BOOL=false', default='false', dest='pairedEnd')
	parser.add_option('-w', '--weighted', help='Specify whether to weight each sample by class distribution', metavar='BOOL=false', default='false', dest='weight')
	parser.add_option('-r', '--normalize_kmer', help='Specify whether or not to normalize kmer counts.  0 = no normalization (useful for assembled data), 1 = normalize by total number of kmers (fast), 2 = normalize by Markov model style.', metavar='INT[0=none,1=by total,2=Markov model]=0', type=int, default=0, dest='normalize')
	parser.add_option('-R', '--early_stop', help='Specify the early stopping rounds', metavar='INT=25', type=int, default=25, dest='earlyStop')
	parser.add_option('-I', '--protein_family', help='Run a conserved protein family/gene model.  Must specify fasta to predict and list of conserved protein families with average lengths, and the number of PLFs to use separated by commas', metavar='FASTA,PLF_List=', default='', dest='runPLF')
	parser.add_option('-v', '--iterative_training', help='Specify whether or not to use interative training with XGBoost.  -1 means no iterative training, otherwise specify batch size to use.', metavar='INT=-1', type=int, default=-1, dest='batchSz')
	parser.add_option('-L', '--alignment', help = 'Run a model using alignment data.  Alignments must be one hot encoded and be in a 2 column tab delimted file of Genome ID, alignment.  No delimiter in alignment', metavar = 'FILE=', default = '', dest='alignmentFile')
	parser.add_option('-u', '--cluster_weight', help = 'Weight each sample by weighting specified in a cluster.  0 = no cluster weight, 1 = weight by size, 2 = weight by SIR cluster dist, 3 = weight by both', metavar='INT=0', default=0, type=int, dest='clustWeight')
	parser.add_option('-U', '--cluster_file', help = 'File to use with cluster weighting (-u|--cluster_weight) option', metavar='FILE=""', default='', dest = 'clustFile')
	parser.add_option('-W', '--two_fold', help = '2 Fold dilution factor values?', metavar = 'BOOL="TRUE"', default = 'True', dest = 'twoFold')

	# parser.add_option('-e', '--cache_kmc', help='Specify whether or not to cache the KMC data for all the genomes in memory between antibiotics', metavar='BOOL=false', default='false', dest='cacheKMC')

	options,args = parser.parse_args()

	options.public = strToBool(options.public)
	options.presence_absence = strToBool(options.presence_absence)
	options.individualModels = strToBool(options.individualModels)
	options.enumerateClasses = strToBool(options.enumerateClasses)
	options.svns = strToBool(options.svns)
	options.noI = strToBool(options.noI)
	options.classify = strToBool(options.classify)
	options.cleanup = strToBool(options.cleanup)
	options.statsOnly = strToBool(options.statsOnly)
	options.pairedEnd = strToBool(options.pairedEnd)
	options.weight = strToBool(options.weight)
	options.header = strToBool(options.header)
	options.twoFold = strToBool(options.twoFold)
	# options.cacheKMC = strToBool(options.cacheKMC)

	return options, parser

def checkOptions(options, parser):
	flag = False

	if options.public:
		if options.fastaDir == '' and options.alignmentFile == '':
			print "fasta directory (-f|--fasta_dir) or alignment file (-L|--alignment) required"
			flag = False
		if options.tabFile == '':
			print "tabular file (-t|--tabular_file) required"
			flag = False
	
		flag = checkDir(options.fastaDir)
		flag = checkFile(options.tabFile)

	if options.tempDir == '':
		tempDir = '.'

	if options.outDir == '':
		outDir = 'model'

	options.fastaDir = cleanDir(options.fastaDir)
	options.tempDir = cleanDir(options.tempDir)
	options.outDir = cleanDir(options.outDir)
	options.kmcDir = cleanDir(options.kmcDir)

	if options.kmcDir != '':
		flag = checkDir(options.kmcDir)
	if options.fastaDir != '':
		flag = checkDir(options.fastaDir)

	if options.statsOnly:
		flag = makeDir(options.outDir, rm = False)
		flag = makeDir(options.tempDir, rm = False)
	else:
		flag = makeDir(options.outDir)
		flag = makeDir(options.tempDir)

	if options.breakpointsFile != '':
		flag = checkFile(options.breakpointsFile)

	if options.totalFolds < options.foldsToRun:
		print "Total folds (-A|--total_folds) must be > folds to run (-a|--folds_to_run)"
		flag = False

	if options.parameters[0] != '{' and options.parameters[-1] != '}':
		print "Properly formatted hash required for model parameters"
		flag = False
	else:
		options.parameters = literal_eval(options.parameters)

	if options.breakpointsFile != '' and options.genus == '':
		print "Warning: Use of breakpoints without a genus is not recommended"
		print "Check to ensure that breakpoints file is correct ONLY for the genus trained on!"

	if flag == False:
		parser.print_help()
		exit(1)

options,parser = getOptions()
checkOptions(options, parser)

print options
os.environ["OMP_NUM_THREADS"] = str(options.threads)

tab = []
gids = {}
labs = {}
if options.public:
	tab, gids, labs = PATRICPublicTabular.parseTabular(options)
else:
	tab, gids, labs = Tabular.loadTab(options)

f = open(options.outDir + 'temp.txt', 'w')
f.write(str(labs) + '\n')
f.close()

# tab = Tabular.shuffleTab(tab, options)

if options.runPLF != '':
	PLF.getPLFFastas(options)

if not options.statsOnly:
	allFeats = {}
	if options.alignmentFile != '':
		allFeats = Alignment.getAllFeats(options)
	else:
		if options.kmcDir == '':
			options.kmcDir = options.tempDir + 'kmc/'
			if os.path.exists('/dev/shm'):
				options.kmcDir = '/dev/shm/kmc' + str(os.getpid()) + '/'
			if not makeDir(options.kmcDir):
				print "Could not make directory:", options.kmcDir
				exit(2)
			KMC.runKMC(options, gids)
		KMC.mergeFastasAndRunKMC(options)
		allFeats = KMC.getAllFeats(options)

	count = 0
	for i in sorted(allFeats):
		allFeats[i] = count
		count += 1

antibioticList = {}
methodsList = {}
for i in tab:
	if i[1] not in antibioticList:
		antibioticList[i[1]] = 0
	if i[2] not in methodsList:
		methodsList[i[2]] = 0
antibioticList = list(sorted(antibioticList))
methodsList = list(sorted(methodsList))

if not options.statsOnly:
	for i in sorted(antibioticList):
		allFeats[i] = len(allFeats)
	for i in sorted(methodsList):
		allFeats[i] = len(allFeats)

tabHsh = Tabular.splitByAntibiotic(tab, options)

if options.stats == '':
	options.stats = []
else:
	options.stats = options.stats.split(',')

for i in tabHsh:
	if options.individualModels and i == 'all':
		continue

	if not options.statsOnly:
		makeDir(options.outDir +  i + '/', rm = False)
		LibSVM.saveLibSVM(options, tabHsh[i], allFeats, i)
		LibSVM.makeAndMergeFolds(options, i)

		weights = None
		if options.weight or options.clustWeight > 0:
			weights = LibSVM.getWeights(options, i)
			LibSVM.writeWeights(options, weights)

		if options.modelType == "XGBoost":
			tab = []
			libsvm = []
			tab, libsvm = LibSVM.loadLibSVM(options, i, weight=weights)

			evals_result = {}
			TrainXGBoost.train(options, tab, labs, libsvm, i)
		else:
			TrainSciKit.train(options, tab, labs, allFeats, weights, i)

	if len(options.stats) != 0:
		for j in options.stats:
			if j == 'rawAcc':
				ModelMetrics.getRawAcc(options, i)
			elif j == 'w1Acc':
				ModelMetrics.getW1Acc(options, i)
			elif j == 'VMEME':
				ModelMetrics.getVMEME(options, i)
			elif j == 'classReport':
				ModelMetrics.classificationReport(options, labs, i)
			elif j == 'confMatrix':
				ModelMetrics.confusionMatrix(options, labs, i)
			elif j == 'f1':
				ModelMetrics.f1(options, labs, i)
			elif j == 'r2':
				ModelMetrics.computeR2(options, i)
			elif j == 'AMRcls':
				ModelMetrics.getRawAcc(options,i)
				ModelMetrics.getVMEME(options, i)
				# ModelMetrics.classificationReport(options, labs, i)
				ModelMetrics.confusionMatrix(options, labs, i)
				ModelMetrics.f1(options, labs, i)
			elif j == 'AMRreg':
				ModelMetrics.getW1Acc(options, i)
				ModelMetrics.getVMEME(options, i)
			elif j == 'cls':
				ModelMetrics.classificationReport(options, labs, i)
				ModelMetrics.confusionMatrix(options, labs, i)
				# ModelMetrics.f1(options, labs, i)
			else:
				print j, "Metric computation is unsupported"

if options.individualModels:
	for i in range(0,options.foldsToRun):
		fList = glob.glob(options.outDir + '*/model.' + str(i) + '.pred')
		fout = open(options.tempDir + 'model.' + str(i) + '.pred', 'w')
		for j in fList:
			f = open(j)

			for k in f:
				fout.write(k)

			f.close()
		fout.close()

		fList = glob.glob(options.outDir + '*/model.' + str(i) + '.true')
		fout = open(options.tempDir + 'model.' + str(i) + '.true', 'w') 
		for j in fList:
			f = open(j)

			for k in f:
				fout.write(k)

			f.close()
		fout.close()

	os.system('mkdir ' + options.outDir + 'all')
	os.system('cp ' + options.tempDir + 'model.*.true ' + options.outDir + 'all/')
	os.system('cp ' + options.tempDir + 'model.*.pred ' + options.outDir + 'all/')

	if len(options.stats) != 0:
		for j in options.stats:
			if j == 'rawAcc':
				ModelMetrics.getRawAcc(options, i)
			elif j == 'w1Acc':
				ModelMetrics.getW1Acc(options, i)
			elif j == 'VMEME':
				ModelMetrics.getVMEME(options, i)
			elif j == 'classReport':
				ModelMetrics.classificationReport(options, labs, i)
			elif j == 'confMatrix':
				ModelMetrics.confusionMatrix(options, labs, i)
			elif j == 'f1':
				ModelMetrics.f1(options, labs, i)
			elif j == 'r2':
				ModelMetrics.computeR2(options, i)
			elif j == 'AMRcls':
				ModelMetrics.getRawAcc(options,i)
				ModelMetrics.getVMEME(options, i)
				# ModelMetrics.classificationReport(options, labs, i)
				ModelMetrics.confusionMatrix(options, labs, i)
				ModelMetrics.f1(options, labs, i)
			elif j == 'AMRreg':
				ModelMetrics.getW1Acc(options, i)
				ModelMetrics.getVMEME(options, i)
			elif j == 'cls':
				ModelMetrics.classificationReport(options, labs, i)
				ModelMetrics.confusionMatrix(options, labs, i)
				# ModelMetrics.f1(options, labs, i)
			else:
				print j, "Metric computation is unsupported"

if options.runPLF != '':
	sDir = os.path.dirname(os.path.realpath(sys.argv[0]))
	cmd = sDir + 'predict.py -f ' + options.tempDir + 'temp.filt.fasta' + ' -m ' + options.outDir + ' -T ' + options.tempDir + ' > ' + options.outDir + 'plf.fasta.preds.tab'
	os.system(cmd)

if not options.statsOnly:
	ModelCleanup.moveFiles(options, gids, allFeats, i, labs, tabHsh)

