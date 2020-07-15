#!/usr/bin/env python

'''
python predict2.py [-f fasta_file] [-m model_dir] <-s species> <-T temp_dir> <-t threads>
'''

from sys import stderr
from optparse import OptionParser
import os
# import KMC
# import Tabular
# import PATRICPublicTabular
import xgboost as xgb
# import Matrix
# import LibSVM
# import TrainXGBoost
from ast import literal_eval
# import ModelCleanup
# import ModelMetrics
# import TrainSciKit
import time
import glob

stderr.write('PID: ' + str(os.getpid()) + '\n')

def getOptions():
	parser = OptionParser()

	parser.add_option('-f', '--fasta', help='Specify fasta file to make predictions for', metavar='FILE', default='', dest='fastaFile')
	parser.add_option('-T', '--temp_dir', help='Specify temporary directory to use', metavar='DIR', default='temp', dest='tempDir')
	parser.add_option('-m', '--model_dir', help='Specify the model directory created from buildModel.py', metavar='MODEL_DIR', default='', dest='modelDir')
	parser.add_option('-t', '--threads', help='Specify the number of threads', metavar='INT=1', type='int', default=1, dest='threads')
	parser.add_option('-i', '--ignore_blank', help='Specify whether or not to ignore blanks in attr order', metavar='BOOL=True', default='True', dest='ignoreBlank')
	parser.add_option('-s', '--species', help='Specify the species being used.  Required for multispecies models', metavar='SPECIES=""', default='', dest='spc')

	options,args = parser.parse_args()

	flag = True

	if options.fastaFile == '':
		flag = False
		print "Fasta file required"
	if options.modelDir == '':
		flag = False
		print "Model directory required"

	if not flag:
		parser.print_help()
		exit(1)

	if options.modelDir[-1] != '/':
		options.modelDir += '/'
	if options.tempDir[-1] != '/':
		options.tempDir += '/'

	if not os.path.isdir(options.tempDir):
		os.mkdir(options.tempDir)

	if options.ignoreBlank[0] == 'T':
		options.ignoreBlank = True
	else:
		options.ignoreBlank = False

	options.spc.replace(' ', '_', len(options.spc))

	return options, parser

def getModelOptions(options):
	f = open(options.modelDir + 'model.params')

	params = f.readline().strip('\n')

	f.close()

	badChars = ';=+*()'
	flag = True
	for i in badChars:
		if i in params:
			flag = False
			print i
			break
	if not flag:
		print "Ill-formated model params"
		exit(1)

	params = literal_eval(params)

	return params

def runKMC(modOpts, options, fName):
	if not os.path.isdir(options.tempDir + 'kmc/'):
		os.mkdir(options.tempDir + 'kmc/')

	cmdArr = []
	tDir = options.tempDir
	# os.mkdir(tDir)
	if 'pairedEnd' in modOpts:
		if modOpts['pairedEnd']:
			cmdArr = ['kmc.sh', str(modOpts['kmerSize']), '@' + fName, options.tempDir + 'kmc/' + os.path.basename(fName), tDir]
		else:
			cmdArr = ['kmc.sh', str(modOpts['kmerSize']), fName, options.tempDir + 'kmc/' + os.path.basename(fName), tDir]
	else:
		cmdArr = ['kmc.sh', str(modOpts['kmerSize']), fName, options.tempDir + 'kmc/' + os.path.basename(fName), tDir]

	cmd = ' '.join(cmdArr) + '> /dev/null'

	os.system(cmd)

def normalizeKmersByTot(kmerHsh):
	tot = 0

	for i in kmerHsh:
		tot += kmerHsh[i]

	tot = float(tot)

	for i in kmerHsh:
		kmerHsh[i] /= tot

def normalizeKmersByMarkov(kmerHsh):
	subHsh = {}

	for i in kmerHsh:
		sub = i[:len(i)/2]
		if sub not in subHsh:
			subHsh[sub] = 0
		subHsh[sub] += kmerHsh[i]

		sub = i[len(i)/2:]
		if sub not in subHsh:
			subHsh[sub] = 0
		subHsh[sub] += kmerHsh[i]

	for i in subHsh:
		subHsh[i] = float(subHsh[i])

	for i in kmerHsh:
		sub = i[:len(i)/2]
		kmerHsh[i] /= subHsh[sub]

def readKMC(modOpts, options, fName, attrOrder):
	kmerHsh = {}
	f = open(options.tempDir + 'kmc/' + os.path.basename(fName) + '.' + str(modOpts['kmerSize']) + '.kmrs')

	for i in f:
		i = i.strip('\n').split('\t')

		if modOpts['presence_absence']:
			kmerHsh[i[0]] = 1
		else:
			kmerHsh[i[0]] = int(i[1])

	f.close()

	if 'normalize' in modOpts:
		if not modOpts['presence_absence']:
			if modOpts['normalize'] == 1:
				normalizeKmersByTot(kmerHsh)
			elif modOpts['normalize'] == 2:
				normalizeKmersByMarkov(kmerHsh)

	kmerStr = ''
	for i in kmerHsh:
		if i not in attrOrder:
			continue
			
		kmerStr += str(attrOrder[i]) + ':' + str(kmerHsh[i]) + ' '

	# return kmerHsh
	return kmerStr

def isKmer(modOpts, st):
	if len(st) != modOpts['kmerSize']:
		return False

	bases = set(['a', 'c', 't', 'g'])

	for i in st:
		if i.lower() not in bases:
			return False

	return True

def parseAttrOrder(modOpts, options):
	f = open(options.modelDir + 'model.attrOrder')

	attrOrder = {}
	antibioticList = []
	for i in f:
		i = i.strip('\n').split('\t')

		if options.ignoreBlank:
			if i[0] == '':
				continue

		attrOrder[i[0]] = int(i[1])

		if not isKmer(modOpts, i[0]):
			antibioticList.append(i[0])

	f.close()

	return attrOrder, antibioticList

def getAntibioticStr(antibiotic, attrOrder):
	abStr = ''
	for i in antibiotic:
		abStr += str(attrOrder[i]) + ':1 '

	return abStr

def makeLibSVM(tab, attrOrder, modOpts, options):
	currFile = ''
	kmerStr = ''

	f = open(options.tempDir + 'sub.libsvm', 'w')
	for i in tab:
		if i[0] != currFile:
			kmerStr = readKMC(modOpts, options, i[0], attrOrder)
			currFile = i[0]

		aStr = getAntibioticStr(i[1:], attrOrder)
		f.write(kmerStr + aStr + '\n')

	f.close()

def predict(modOpts, options):
	model = None
	
	mFile = ""
	if os.path.exists(options.modelDir + 'all/model.0pkl'):
		mFile = options.modelDir + 'all/model.0pkl'
	else:
		mFile = options.modelDir + 'all/model.0.pkl'

	model = xgb.Booster(model_file = mFile)
	model.set_param('nthread', options.threads)

	dm = xgb.DMatrix(options.tempDir + 'sub.libsvm', silent = True)
	pred = model.predict(dm)

	return pred

def getLabelMap(options):
	f = open(options.modelDir + 'model.labels.map')

	labelMap = {}
	for i in f:
		i = i.strip('\n').split('\t')
		# stderr.write(str(i) + '\n')
		try:
			labelMap[int(i[1])] = i[0]
		except:
			labelMap[float(i[1])] = i[0]

	f.close()

	return labelMap

def getSpcAbCombos(options):
	f = open(options.modelDir + 'spc_ab_std.tab')

	antibioticList = []
	for i in f:
		i = i.strip('\n').split('\t')
		if options.spc == i[0]:
			antibioticList.append(i[1:])

	f.close()

	return antibioticList

options, parser = getOptions()
modOpts = getModelOptions(options)

attrOrder, antibioticList = parseAttrOrder(modOpts, options)

if options.spc != '' and os.path.exists(options.modelDir + 'spc_ab_std.tab'):
	antibioticList = getSpcAbCombos(options)

labelMap = getLabelMap(options)

fileList = []
if os.path.isdir(options.fastaFile):
	if options.fastaFile[-1] != '/':
		options.fastaFile += '/'
	fileList = glob.glob(options.fastaFile + '*.fasta')
else:
	fileList = [options.fastaFile]

tab = []
count = 1
inc = len(fileList) / 50
stderr.write("Making predictions...\n\t")
for i in fileList:
	# if count % inc == 0:
	# 	stderr.write('=')

	runKMC(modOpts, options, i)
	for j in antibioticList:
		if isinstance(j, list):
			tab.append([i] + j)
		else:
			tab.append([i, j])

	if count % (options.threads*2) == 0:
		makeLibSVM(tab, attrOrder, modOpts, options)
		pred = predict(modOpts, options)

		for j in range(0,len(tab)):
			tab[j].append(str(pred[j]))
			pred[j] = round(pred[j])
			if int(pred[j]) in labelMap:
				tab[j].append(labelMap[int(pred[j])])
			print '\t'.join(tab[j])

		tab = []

		# break

	count += 1

if len(tab) > 0:
	makeLibSVM(tab, attrOrder, modOpts, options)
	pred = predict(modOpts, options)

	for j in range(0,len(tab)):
		tab[j].append(str(pred[j]))
		pred[j] = round(pred[j])
		if int(pred[j]) in labelMap:
			tab[j].append(labelMap[int(pred[j])])
		print '\t'.join(tab[j])

stderr.write('\nDone\n')

