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
from copy import deepcopy

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

def cleanTemp(options):
	os.system('rm -r ' + options.tempDir + '*')

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

def getAbStdComboFromTab(fName):
	f = open(fName)

	combos = {}
	for i in f:
		i = i.strip('\n').split('\t')
		combo = i[1:-1]
		combo = ':'.join(combo)
		if combo not in combos:
			combos[combo] = 0

	f.close()

	return combos

def mrgHsh(to, fr):
	for i in fr:
		if i not in to:
			to[i] = 0

def getAbStdCombos(options):
	fLst = glob.glob(options.modelDir + 'all/model.*.true')

	aCombos = {}
	for i in fLst:
		combos = getAbStdComboFromTab(i)
		mrgHsh(aCombos, combos)

	return aCombos

def runKMC(options, params):
	k = params['kmerSize']
	oFile = options.tempDir + os.path.basename(options.fastaFile)

	cmdArr = ['kmc.sh', str(k), options.fastaFile, oFile, options.tempDir, '> /dev/null']
	cmd = ' '.join(cmdArr)

	os.system(cmd)

	return oFile + '.' + str(k) + '.kmrs'

def getAttrOrder(options):
	f = open(options.modelDir + 'model.attrOrder')

	attrOrder = {}
	for i in f:
		i = i.strip('\n').split('\t')
		attrOrder[i[0]] = int(i[1])

	f.close()

	return attrOrder

def getKmerArr(attrOrder, kmcFile):
	arr = [0] * len(attrOrder)

	f = open(kmcFile)

	for i in f:
		i = i.strip('\n').split('\t')
		if i[0] not in attrOrder:
			stderr.write(kmcFile + '\t' + i[0] + '\n')
			continue
		ind = attrOrder[i[0]]
		arr[ind] = int(i[1])

	f.close()

	return arr

def makeMatrix(arr, attrOrder, aCombos):
	mat = []
	cOrd = []
	for i in sorted(aCombos):
		cOrd.append(i)
		i = i.split(':')
		arrCopy = deepcopy(arr)
		for j in i:
			ind = attrOrder[j]
			if ind >= len(arrCopy):
				continue
			arrCopy[ind] = 1

		mat.append(arrCopy)

	return mat, cOrd

def printLibSVM(mat, options):
	mStr = ''
	for i in mat:
		lStr = '0 '
		for j in range(0,len(i)):
			if i[j] != 0:
				lStr += ' ' + str(j) + ':' + str(i[j])
		mStr += lStr[:-1] + '\n'

	fout = options.tempDir + 'mat.libsvm'
	f = open(fout, 'w')
	f.write(mStr)
	f.close()

	return fout

def getLabMap(options):
	f = open(options.modelDir + 'model.labels.map')

	labHsh = {}
	for i in f:
		i = i.strip('\n').split('\t')
		labHsh[float(i[1])] = i[0]

	f.close()

	return labHsh

def printPreds(preds, options, aCombos):
	labHsh = getLabMap(options)
	fName = os.path.basename(options.fastaFile)

	for i in range(0,len(preds)):
		pred = str(preds[i])
		lab = ''
		if preds[i] not in labHsh:
			lab = str(round(preds[i], 0))
		else:
			lab = labHsh[preds[i]]
		combo = aCombos[i]

		arr = [fName, combo, pred, lab]

		print '\t'.join(arr)

def main():
	options,parser = getOptions()
	cleanTemp(options)
	params = getModelOptions(options)

	aCombos = getAbStdCombos(options)
	attrOrder = getAttrOrder(options)

	kmcFile = runKMC(options, params)
	arr = getKmerArr(attrOrder, kmcFile)

	mat,cOrd = makeMatrix(arr, attrOrder, aCombos)
	lFile = printLibSVM(mat, options)

	dMat = xgb.DMatrix(lFile, silent = True)
	mod = xgb.Booster(model_file = options.modelDir + 'all/model.0pkl')

	preds = mod.predict(dMat)

	printPreds(preds, options, cOrd)

if __name__ == '__main__':
	main()
