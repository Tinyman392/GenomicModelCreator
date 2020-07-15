'''
python rawOutTable.py [model directory] [breakpoints,genus] [antibiotic]
'''

from sys import argv
from glob import glob
import os
from math import log

if argv[1][-1] != '/':
	argv[1] += '/'

def parsePredFile(fName):
	f = open(fName)

	pred = []
	for i in f:
		i = i.strip()
		pred.append(round(float(i)))

	f.close()

	return pred

def parseTrueFile(fName):
	f = open(fName)

	true = []
	for i in f:
		i = i.strip().split('\t')
		gid = i[0]
		antibiotic = i[1].split(':')[0]
		label = round(float(i[2]))

		true.append([gid, antibiotic, label])

	f.close()

	return true

def parseModelDir(modelDir):
	trueFs = glob(modelDir + '*.true')
	predFs = glob(modelDir + '*.pred')

	true = []
	pred = []

	for  i in range(0,len(trueFs)):
		true = true + parseTrueFile(trueFs[i])
		pred = pred + parsePredFile(predFs[i])

	return true, pred

def getModelDir():
	modelDir = argv[1]
	if len(argv) > 3:
		if os.path.isdir(modelDir + argv[3] + '/'):
			modelDir = modelDir + argv[3] + '/'

	return modelDir

def parseBreaksFile(breaksFile, genus):
	f = open(breaksFile)
	f.readline()

	breaks = {}
	for i in f:
		if len(i) == 0 or i[0] == '#':
			continue
		i = i.strip().split('\t')

		antibiotic = ''
		sBreak = 0
		rBreak = 0

		if len(i) < 5:
			antibiotic = i[0]
			sBreak = round(log(float(i[1]), 2))
			rBreak = round(log(float(i[2]),2))

			breaks[antibiotic] = [sBreak, rBreak]
		else:
			if i[0] != genus:
				continue
			antibiotic = i[1]
			sBreak = round(log(float(i[2]), 2))
			rBreak = round(log(float(i[3]),2))

			breaks[antibiotic] = [sBreak, rBreak]

	f.close()

	return breaks

def getBreakpoints(true):
	breaks = {}
	if len(argv) > 2:
		breaksFile = argv[2].split(',')[0]
		genus = argv[2].split(',')[-1]

		if genus == breaksFile:
			genus = ''

		breaks = parseBreaksFile(breaksFile, genus)
	else:
		for i in true:
			antibiotic = i[1]
			if antibiotic not in breaks:
				breaks[antibiotic] = [0, 0]
			if i[-1] > breaks[antibiotic][-1]:
				breaks[antibiotic][-1] = round(log(i[-1],2))
			elif i[-1] < breaks[antibiotic][0]:
				breaks[antibiotic][0] = round(log(i[0], 2))

	return breaks

def merge(true, pred, breaks):
	tab = []
	for i in range(0,len(true)):
		t = true[i][-1]
		a = true[i][1]
		p = pred[i]

		err = '0'

		if t != p:
			if a in breaks:
				if p <= breaks[a][0] and t >= breaks[a][1]:
					err = '2'
				elif p >= breaks[a][1] and t <= breaks[a][0]:
					err = '1'

		arr = true[i] + [p, err]

		tab.append(arr)

	for i in sorted(tab, key = lambda x: x[-1], reverse = True):
		if i[-1] == '0':
			i[-1] = ''
		elif i[-1] == '1':
			i[-1] = 'ME'
		elif i[-1] == '2':
			i[-1] = 'VME'

		for j in range(0,len(i)):
			i[j] = str(i[j])

		print '\t'.join(i)

modelDir = getModelDir()

true, pred = parseModelDir(modelDir)
breaks = getBreakpoints(true)

merge(true, pred, breaks)


