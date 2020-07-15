import KMC
from sys import stderr
import os
import xgboost as xgb
from math import log
from random import randint
import Alignment

def toMIC(micStr):
	mic = 0
	for i in range(0,len(micStr)):
		try:
			mic = float(micStr[i:])
			break
		except:
			continue

	if micStr[0] == '>' and '=' not in micStr:
		mic *= 2
	elif micStr[0] == '<' and '=' not in micStr:
		mic /= 2

	# mic = log(mic,2)

	return int(round(mic,0))

def kmerHshToKmerStr(kmerHsh, allFeats):
	kmerStr = ''
	for i in kmerHsh:
		ind = -1
		try:
			ind = allFeats[i]
		except:
			allFeats[i] = len(allFeats)
			ind = allFeats[i]
		
		kmerStr += ' ' + str(ind) + ':' + str(kmerHsh[i])

	return kmerStr

def listToLibSVMStr(lst):
	st = ''
	for i in range(0,len(lst)):
		st += ' ' + str(i) + ':' + str(lst[i])

	return st

def saveLibSVM(options, tab, allFeats, fLabel='all'):
	if options.alignmentFile == '':
		kmerStr = ''
		currGID = ''

		fLSVM = open(options.tempDir + fLabel + '.libsvm', 'w')
		fTrue = open(options.tempDir + fLabel + '.true', 'w')

		stderr.write("Writing LibSVM...\n\t")
		count = 0
		inc = len(tab) / 50
		for i in tab:
			if count > inc:
				stderr.write('=')
				count = 0
			count += 1

			gid = i[0]
			if currGID != gid:
				currGID = gid
				kmerHsh = KMC.readKMCOut(options.kmcDir + gid + '.fasta.' + str(options.kmerSize) + '.kmrs', options)
				KMC.normalizeKMC(kmerHsh, options)
				kmerStr = kmerHshToKmerStr(kmerHsh, allFeats)

			if kmerStr == '':
				continue

			antibiotic = i[1].split(':')[0]
			antibioticStr = ' ' + str(allFeats[antibiotic]) + ':1'
			label = str(i[2])

			lsvmStr = label + kmerStr + antibioticStr + '\n'
			trueStr = i[0] + '\t' + i[1] + '\t' + str(i[2]) + '\n'

			fLSVM.write(lsvmStr)
			fTrue.write(trueStr)
		stderr.write('\n')

		fLSVM.close()
		fTrue.close()
	else:
		print "Getting alignments"
		alignmentHsh = Alignment.parseAlignmentFile(options)
		for i in alignmentHsh:
			alignmentHsh[i] = listToLibSVMStr(alignmentHsh[i])

		fLSVM = open(options.tempDir + fLabel + '.libsvm', 'w')
		fTrue = open(options.tempDir + fLabel + '.true', 'w')

		stderr.write("Writing LibSVM...\n\t")
		count = 0
		inc = len(tab) / 50
		for i in tab:
			if count > inc:
				stderr.write('=')
				count = 0
			count += 1

			gid = i[0]

			if gid not in alignmentHsh:
				continue

			antibiotic = i[1].split(':')[0]
			antibioticStr = ' ' + str(allFeats[antibiotic]) + ':1'
			label = str(i[2])
			alignment = alignmentHsh[gid]

			lsvmStr = label + alignment + antibioticStr + '\n'
			trueStr = i[0] + '\t' + i[1] + '\t' + str(i[2]) + '\n'

			fLSVM.write(lsvmStr)
			fTrue.write(trueStr)
		stderr.write('\n')

		fLSVM.close()
		fTrue.close()


def catFolds(options, fLabel='all'):
	foutT = open(options.tempDir + fLabel + '.true', 'w')
	foutL = open(options.tempDir + fLabel + '.libsvm', 'w')

	for i in range(0,options.totalFolds):
		finT = open(options.tempDir + fLabel + '.' + str(i) + '.spl.true')
		finL = open(options.tempDir + fLabel + '.' + str(i) + '.spl.libsvm')

		for i in finT:
			foutT.write(i)
		for i in finL:
			foutL.write(i)

def makeAndMergeFolds(options, fLabel='all'):
	fTrues = []
	fLSVMs = []
	for i in range(0,options.totalFolds):
		fTrues.append(open(options.tempDir + fLabel + '.' + str(i) + '.spl.true', 'w'))
		fLSVMs.append(open(options.tempDir + fLabel + '.' + str(i) + '.spl.libsvm', 'w'))

	fin = open(options.tempDir + fLabel + '.true')

	comboCounts = {}
	lineToIndHsh = []
	for i in fin:
		line = i
		i = i.strip('\n').split('\t')
		antibiotic = i[1].split(':')[0]
		label = i[2]

		if antibiotic not in comboCounts:
			comboCounts[antibiotic] = {}
		if label not in comboCounts[antibiotic]:
			comboCounts[antibiotic][label] = randint(0, options.totalFolds - 1)

		ind = comboCounts[antibiotic][label] % options.totalFolds
		lineToIndHsh.append(ind)

		comboCounts[antibiotic][label] += 1

		fTrues[ind].write(line)

	fin.close()

	fin = open(options.tempDir + fLabel + '.libsvm')

	stderr.write("Splitting LibSVM...\n\t")
	count = 0
	count2 = 0
	inc = len(lineToIndHsh) / 50
	for i in fin:
		if count2 > inc:
			stderr.write('=')
			count2 = 0
		count2 += 1

		ind = lineToIndHsh[count]
		fLSVMs[ind].write(i)

		count += 1
	stderr.write('\n')

	fin.close()

	for i in range(0,len(fTrues)):
		fTrues[i].close()
		fLSVMs[i].close()

	stderr.write("Catting tabs and libsvms...")
	catFolds(options, fLabel)
	stderr.write("Done\n")

	# stderr.write('Catting tabs...')
	# cmdArr = ['cat', options.tempDir + fLabel + '.*.spl.true', options.tempDir + fLabel + '.true']
	# cmd = ' '.join(cmdArr)
	# os.system(cmd)
	# stderr.write('Done\n')

	# stderr.write('Catting libsvms...')
	# cmdArr = ['cat', options.tempDir + fLabel + '.*.spl.libsvm >', options.tempDir + fLabel + '.libsvm']
	# cmd = ' '.join(cmdArr)
	# os.system(cmd)
	# stderr.write('Done\n')

	stderr.write('Removing split files...')
	cmdArr = ['rm', options.tempDir + fLabel + '*.spl.*']
	cmd = ' '.join(cmdArr)
	os.system(cmd)
	stderr.write('Done\n')

def loadLibSVM(options, fLabel='all', weight = None):
	f = open(options.tempDir + fLabel + '.true')

	tab = []
	for i in f:
		i = i.strip('\n').split('\t')
		i[2] = float(i[2])
		tab.append(i)

	f.close()

	libsvm = xgb.DMatrix(options.tempDir + fLabel + '.libsvm', weight = weight)

	return tab, libsvm

def getWeights(options, fLabel='all'):
	print "Getting weights"
	f = open(options.tempDir + fLabel + '.true')

	tab = []
	counts = {}
	tots = {}
	for i in f:
		i = i.strip().split('\t')
		if i[1] not in counts:
			counts[i[1]] = {}
			tots[i[1]] = 0
		if i[2] not in counts[i[1]]:
			counts[i[1]][i[2]] = 0

		counts[i[1]][i[2]] += 1
		tots[i[1]] += 1

		tab.append(i)

	f.close()

	for i in counts:
		for j in counts[i]:
			counts[i][j] = 1 - float(counts[i][j]) / tots[i]

	f = open(options.tempDir + fLabel + '.weights', 'w')
	weights = []
	for i in tab:
		w = counts[i[1]][i[2]]
		f.write(str(w) + '\n')
		weights.append(w)
	f.close()

	return weights
