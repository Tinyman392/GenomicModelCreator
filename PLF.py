from glob import glob
from sys import stderr
import os

def parsePLFTab(options):
	f = open(options.runPLF.split(',')[1].strip())

	plfTab = {}
	for i in f:
		i = i.strip().split('\t')
		plfTab[i[0]] = float(i[1])

	f.close()

	return plfTab

def parseFasta(fName):
	f = open(fName)

	currSeq = ''
	currGID = ''
	fastaHsh = {}
	for i in f:
		i = i.strip()
		if len(i) == 0:
			continue
		if i[0] == '>':
			if currSeq != '' and currGID != '':
				if currGID not in fastaHsh:
					fastaHsh[currGID] = []
				fastaHsh[currGID].append(currSeq)

			currSeq = ''
			currGID = i.split('>')[1].split(' ')[0]

		currSeq += i

	if currSeq != '' and currGID != '':
		if currGID not in fastaHsh:
			fastaHsh[currGID] = []
		fastaHsh[currGID].append(currSeq)

	f.close()

	return fasta

def getSumLen(arr):
	tot = 0
	for i in arr:
		tot += len(i)

	return tot

def filterFasta(plfList, fastaHsh):
	delList = []
	for i in fastaHsh:
		if i not in plfList:
			delList.append(i)

	for i in delList:
		del fastaHsh[i]

def getTopPLFs(options, plfTab, fastaHsh):
	plfList = []
	for i in fastaHsh:
		if i in plfTab:
			tLen = getSumLen(fastaHsh[i])
			if tLen >= plfTab[i] / 2 and tLen <= plfTab[i] * 2:
				plfList.append(i)

	nPLF = options.runPLF.split(',')[3]
	if nPLF.lower() != 'max':
		nPLF = int(nPLF)
		plfList.sort(key = lambda x: plfTab[x], reverse = True)[:nPLF]

		if len(plfList) < nPLF:
			stderr.write("Not enough PLFs to build model...\n")
			exit(1)

	plfList = set(plfList)

	filterFasta(plfList, fastaHsh)

	return plfList

def writeFasta(fastaHsh, fout):
	f = open(fout, 'w')

	for i in fastaHsh:
		for j in fastaHsh[i]:
			f.write('>' + i + '\n')
			f.write(j + '\n')

	f.close()

def filterFastas(options, plfList):
	fList = glob(options.fastaDir + '*.fasta')
	oDir = options.tempDir + 'fasta/'

	for i in fList:
		fastaHsh = parseFasta(i)
		filterFasta(plfList, fastaHsh)
		base = os.path.basename(i)
		writeFasta(fastaHsh, oDir + base)

	options.fastaDir = oDir

def writePLFList(options, plfList):
	f = open(options.outDir + 'plfs.used.list', 'w')

	for i in plfList:
		f.write(i + '\n')

	f.close()

def getPLFFastas(options):
	plfTab = parsePLFTab(options)
	fastaHsh = parseFasta(options.runPLF.split(',')[0].strip())
	plfList = getTopPLFs(options, plfTab, fastaHsh)
	filterFastas(options, plfList)

	writeFasta(fastaHsh, options.tempDir + 'temp.filt.fasta')
	writePLFList(options, plfList)

