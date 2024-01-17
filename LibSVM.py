import KMC
from sys import stderr
import os
import xgboost as xgb
from math import log, isnan
from random import randint
import Alignment
import TabFeat

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

def parseDrugDesc(fName, nKmers):
	if fName == '':
		return {}
	f = open(fName)

	descHsh = {}
	for i in f:
		i = i.strip('\n').split('\t')
		ab1 = i[1]
		ab2 = i[0]
		arr = i[2].split(',')
		descStr = ''

		for j in range(0,len(arr)):
			if 'nan' in arr[j]:
				continue
			descStr += ' ' + str(j+nKmers) + ':' + arr[j]

		descHsh[ab1] = descStr
		descHsh[ab2] = descStr

	f.close()

	return descHsh

def saveLibSVM(options, tab, allFeats, fLabel='all'):
	descHsh = parseDrugDesc(options.drugDescFile, len(allFeats))
	if options.alignmentFile == '' and options.tabFeat == '':
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

				if not os.path.exists(options.kmcDir + gid + '.fasta.' + str(options.kmerSize) + '.kmrs'):
					kmerStr = ''
					stderr.write("KMC doesn't exist for: " + gid + '\n')
					continue

				kmerHsh = KMC.readKMCOut(options.kmcDir + gid + '.fasta.' + str(options.kmerSize) + '.kmrs', options)
				KMC.normalizeKMC(kmerHsh, options)
				kmerStr = kmerHshToKmerStr(kmerHsh, allFeats)

			if kmerStr == '':
				continue

			antibiotic = i[1]
			if len(descHsh) == 0:
				antibioticStr = ' ' + str(allFeats[antibiotic]) + ':1'
			else:
				if antibiotic in descHsh:
					antibioticStr = descHsh[antibiotic]
				else:
					continue
			method = i[2]
			methodStr = ' ' + str(allFeats[method]) + ':1'
			label = str(i[3])

			lsvmStr = label + kmerStr + antibioticStr + methodStr + '\n'
			trueStr = i[0] + '\t' + i[1] + '\t' + i[2] + '\t' + str(i[3]) + '\n'

			fLSVM.write(lsvmStr)
			fTrue.write(trueStr)
		stderr.write('\n')

		fLSVM.close()
		fTrue.close()
	elif options.alignmentFile == '':
		stderr.write("Parsing tabular features\n")
		tabFeatHsh = TabFeat.parseFeats(options)

		stderr.write('Converting Array to LibSVM Format...\n\t')
		count = 0
		inc = len(tabFeatHsh) / 50
		for i in tabFeatHsh:
			if count > inc:
				stderr.write('=')
				count = 0
			count += 1

			tabFeatHsh[i] = listToLibSVMStr(tabFeatHsh[i])
		stderr.write('\n')

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

			if gid not in tabFeatHsh:
				continue

			antibiotic = i[1]

			if len(descHsh) == 0:
				antibioticStr = ' ' + str(allFeats[antibiotic]) + ':1'
			else:
				if antibiotic in descHsh:
					antibioticStr = descHsh[antibiotic]
				else:
					continue

			method = i[2]
			methodStr = ' ' + str(allFeats[method]) + ':1'

			label = str(i[3])

			feat = tabFeatHsh[gid]

			lsvmStr = label + feat + antibioticStr + methodStr + '\n'
			trueStr = i[0] + '\t' + i[1] + '\t' + i[2] + '\t' + str(i[3]) + '\n'

			fLSVM.write(lsvmStr)
			fTrue.write(trueStr)
		stderr.write('\n')

		fLSVM.close()
		fTrue.close()
	else:
		print "Getting alignments"
		alignmentHsh = Alignment.parseFeats(options)
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
				# stderr.write(gid + '\n')
				continue

			antibiotic = i[1]
			if len(descHsh) == 0:
				antibioticStr = ' ' + str(allFeats[antibiotic]) + ':1'
			else:
				if antibiotic in descHsh:
					antibioticStr = descHsh[antibiotic]
				else:
					continue
			method = i[2]
			methodStr = ' ' + str(allFeats[method]) + ':1'
			label = str(i[3])

			alignment = alignmentHsh[gid]

			lsvmStr = label + alignment + antibioticStr + methodStr + '\n'
			trueStr = i[0] + '\t' + i[1] + '\t' + i[2] + '\t' + str(i[3]) + '\n'

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

	foutT.close()
	foutL.close()

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
		label = i[3]

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
		i[3] = float(i[3])
		tab.append(i)

	f.close()

	libsvm = xgb.DMatrix(options.tempDir + fLabel + '.libsvm', weight = weight)

	return tab, libsvm

def getWeights(options, fLabel='all'):
	print "Getting weights"
	weights = []

	f = open(options.tempDir + fLabel + '.true')

	tab = []
	tabHsh = {}
	counts = {}
	tots = {}
	for i in f:
		i = i.strip().split('\t')
		if i[1] not in counts:
			counts[i[1]] = {}
			tots[i[1]] = 0
		if i[3] not in counts[i[1]]:
			counts[i[1]][i[3]] = 0

		if i[0] not in tabHsh:
			tabHsh[i[0]] = {}
		tabHsh[i[0]][i[1]] = i[3]

		counts[i[1]][i[3]] += 1
		tots[i[1]] += 1

		tab.append(i)

	f.close()

	if options.clustWeight == 0 and options.weight == True:
		for i in counts:
			for j in counts[i]:
				counts[i][j] = 1 - float(counts[i][j]) / tots[i]

		f = open(options.tempDir + fLabel + '.weights', 'w')
		
		for i in tab:
			w = counts[i[1]][i[3]]
			f.write(str(w) + '\n')
			weights.append(w)
		f.close()
	else:
		f = open(options.clustFile)

		genClustHsh = {}
		clustGenHsh = {}
		for i in f:
			i = i.strip().split('\t')

			genClustHsh[i[0]] = i[1]

			if i[1] not in clustGenHsh:
				clustGenHsh[i[1]] = []
			clustGenHsh[i[1]].append(i[0])

		f.close()		

		weights1 = []
		if options.clustWeight == 1 or options.clustWeight == 3:
			for i in tab:
				clust = genClustHsh[i[0]]
				sz = len(clustGenHsh[clust])
				weights1.append(1 - float(sz) / len(tab))

		weights2 = []
		if options.clustWeight == 2 or options.clustWeight == 3:
			clustAntibioLabHsh = {}
			for i in tab:
				clust = genClustHsh[i[0]]
				if clust not in clustAntibioLabHsh:
					clustAntibioLabHsh[clust] = {}
				if i[1] not in clustAntibioLabHsh[clust]:
					clustAntibioLabHsh[clust][i[1]] = {}
				if i[3] not in clustAntibioLabHsh[clust][i[1]]:
					clustAntibioLabHsh[clust][i[1]][i[3]] = 0
				clustAntibioLabHsh[clust][i[1]][i[3]] += 1

			for i in clustAntibioLabHsh:
				for j in clustAntibioLabHsh[i]:
					s = 0
					for k in clustAntibioLabHsh[i][j]:
						s += clustAntibioLabHsh[i][j][k]
					for k in clustAntibioLabHsh[i][j]:
						if s == 0:
							clustAntibioLabHsh[i][j][k] = 1
						else:
							clustAntibioLabHsh[i][j][k] =  1 - float(clustAntibioLabHsh[i][j][k]) / s

			for i in tab:
				clust = genClustHsh[i[0]]
				weight = clustAntibioLabHsh[clust][i[1]][i[3]]
				weights2.append(weight)

		if options.clustWeight == 1:
			weights = weights1
		elif options.clustWeight == 2:
			weights = weights2
		elif options.clustWeight == 3:
			for i in range(0,len(weights1)):
				weights.append(weights1[i] * weights2[i])
		else:
			weights = None

		if weights is not None:
			for i in range(0,len(weights)):
				if weights[i] == 0:
					weights[i] = 1

	return weights

def writeWeights(options, weights):
	f = open(options.outDir + 'weights.list', 'w')

	if weights is not None:
		for i in weights:
			f.write(str(i) + '\n')

	f.close()
