#!/usr/bin/env python

'''
python predict.py [-f fasta_file] [-m model_dir] <-T temp_dir> <-t threads>
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
	tDir = options.tempDir[:-1] + '_' + str(os.getpid()) + '/'
	os.mkdir(tDir)
	if 'pairedEnd' in modOpts:
		if modOpts['pairedEnd']:
			cmdArr = ['kmc.sh', str(modOpts['kmerSize']), '@' + fName, options.tempDir + 'kmc/' + os.path.basename(fName), tDir]
		else:
			cmdArr = ['kmc.sh', str(modOpts['kmerSize']), fName, options.tempDir + 'kmc/' + os.path.basename(fName), tDir]
	else:
		cmdArr = ['kmc.sh', str(modOpts['kmerSize']), fName, options.tempDir + 'kmc/' + os.path.basename(fName), tDir]

	cmd = ' '.join(cmdArr) + '> /dev/null'

	os.system(cmd)

	os.system('rm -rf ' + tDir)

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

def readKMC(modOpts, options, fName):
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

	return kmerHsh

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
		attrOrder[i[0]] = int(i[1])

		if not isKmer(modOpts, i[0]):
			antibioticList.append(i[0])

	f.close()

	return attrOrder, antibioticList

def kmerHshToArr(kmerHsh, attrOrder):
	arr = [0] * len(attrOrder)
	for i in kmerHsh:
		if i not in attrOrder:
			continue

		arr[attrOrder[i]] = kmerHsh[i]

	return arr

def predict(kmerArr, options, antibiotic='all'):
	model = None
	if os.path.isdir(options.modelDir + antibiotic):
		mFile = ""
		if os.path.exists(options.modelDir + antibiotic + '/model.0pkl'):
			mFile = options.modelDir + antibiotic + '/model.0pkl'
		else:
			mFile = options.modelDir + antibiotic + '/model.0.pkl'
		
		model = xgb.Booster(model_file = mFile)
	else:
		mFile = ""
		if os.path.exists(options.modelDir + 'all/model.0pkl'):
			mFile = options.modelDir + 'all/model.0pkl'
		else:
			mFile = options.modelDir + 'all/model.0.pkl'

		model = xgb.Booster(model_file = mFile)
	model.set_param('nthread', 1)

	fName = options.tempDir + 'temp.child' + str(os.getpid()) + '.libsvm'
	dm = xgb.DMatrix(fName, silent = True)
	dm = xgb.DMatrix(kmerArr, missing = 0)

	pred = model.predict(dm)

	return pred[0]

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

def toKmerStr(kmerHsh, attrOrder):
	kmerStr = '0 '
	for i in kmerHsh:
		if i not in attrOrder:
			continue

		if kmerHsh[i] == -1:
			continue

		ind = attrOrder[i]
		kmerStr += str(ind) + ':' + str(kmerHsh[i]) + ' '

	return kmerStr

def toLibSVM(options, kmerStr, antibiotic, attrOrder):
	f = open(options.tempDir + 'temp.child' + str(os.getpid()) + '.libsvm', 'w')

	f.write(kmerStr + str(attrOrder[antibiotic]) + ':1\n')

	f.close()

def mergeChildOut(options, pidLst):
	# fList = glob.glob(options.tempDir + 'temp.child.out.*')
	fList = []
	for i in pidLst:
		fList.append(options.tempDir + 'temp.child.out.' + str(i))

	for i in fList:
		f = open(i)
		for j in f:
			print j,

		f.close()

	os.system('rm ' + options.tempDir + 'temp.child.out.*')
	os.system('rm ' + options.tempDir + 'kmc/*')
	os.system('rm ' + options.tempDir + '*.libsvm')

	os._exit(0)

options, parser = getOptions()
modOpts = getModelOptions(options)

attrOrder, antibioticList = parseAttrOrder(modOpts, options)
labelMap = getLabelMap(options)

if os.path.isdir(options.fastaFile):
	if options.fastaFile[-1] != '/':
		options.fastaFile += '/'
	fileList = glob.glob(options.fastaFile + '*.fasta')
else:
	fileList = [options.fastaFile]

count = 0
oThreads = 0
start = time.time()
pidLst = []
pbCount = 0
inc = len(fileList) / 50
stderr.write('Making predictions\n\t')
for j in fileList:
	if pbCount >= inc:
		stderr.write('=')
		pbCount = 0
	pbCount += 1

	if oThreads > options.threads - 1:
		for i in pidLst:
			os.waitpid(i, 0)

		mergeChildOut(options, pidLst)
	
		pidLst = []
		oThreads = 0

	pid = os.fork()
	if pid == 0:
		f = open(options.tempDir + 'temp.child.out.' + str(os.getpid()), 'w')
		runKMC(modOpts, options, j)
		kmerHsh = readKMC(modOpts, options, j)
		kmerArr = kmerHshToArr(kmerHsh, attrOrder)
		kmerStr = toKmerStr(kmerHsh, attrOrder)

		for i in antibioticList:
			if i == '':
				continue

			toLibSVM(options, kmerStr, i, attrOrder)

			ind = attrOrder[i]
			kmerArr[ind] = 1
			
			pred = predict(kmerArr, options, i)

			label = pred
			if pred in labelMap:
				label = labelMap[pred]
			else:
				label = str(pred)

			kmerArr[ind] = 0

			f.write('\t'.join([j, i, label]) + '\n')

		f.close()

		os._exit(0)

	pidLst.append(pid)

	oThreads += 1

	# if count > 64:
	# 	break
	# count += 1

for i in pidLst:
	os.waitpid(i,0)

mergeChildOut(options, pidLst)

stderr.write('\n')

end = time.time()

stderr.write(str(end - start) + '\n')

