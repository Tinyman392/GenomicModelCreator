from glob import glob
import os
import re
from sys import stderr

def runKMC(options, gids):
	dList = glob(options.kmcDir + '/*')
	for i in dList:
		os.remove(i)

	fList = []
	for i in gids:
		fList.append(i + '.fasta')

	for i in fList:
		fName = options.fastaDir + i

		cmdArr = []
		if options.pairedEnd:
			cmdArr = ['kmc.sh', str(options.kmerSize), '@' + fName, options.kmcDir + i, options.kmcDir, "> /dev/null"]
		else:
			cmdArr = ['kmc.sh', str(options.kmerSize), fName, options.kmcDir + i, options.kmcDir, "> /dev/null"]
		cmd = ' '.join(cmdArr)

		os.system(cmd)
		os.system('rm ' + str(options.kmcDir) + i + '.kmc_*')

def readKMCOut(fName, options):
	f = open(fName)

	kmerHsh = {}
	for i in f:
		i = i.strip().split('\t')
		if options.presence_absence:
			kmerHsh[i[0]] = 1
		else:
			kmerHsh[i[0]] = int(i[1])

	f.close()

	return kmerHsh

def readKMC(options):
	fList = glob(options.kmcDir + '*.' + str(options.kmerSize) + '.kmrs')

	kmerHsh = {}
	count = 0
	inc = len(fList) / 50
	stderr.write("Reading KMC files...\n\t")
	for i in fList:
		if count > inc:
			stderr.write('=')
			count = 0
		count += 1

		gid = re.sub(r'\.fasta\.[0-9]*\.kmrs', '', os.path.basename(i))
		kmerHsh[gid] = readKMCOut(i, options)
	stderr.write('\n')

	return kmerHsh

def mergeFastasAndRunKMC(options):
	cmdArr = ['cat', str(options.fastaDir) + '*.fasta > ' + str(options.tempDir + 'allFasta.fasta')]
	cmd = ' '.join(cmdArr)
	os.system(cmd)

	if not options.pairedEnd:
		cmdArr = ['kmc.sh', str(options.kmerSize), options.tempDir + 'allFasta.fasta', options.tempDir + 'allFasta.fasta', options.tempDir]
	else:
		cmdArr = ['kmc.sh', str(options.kmerSize), '@' + options.tempDir + 'allFasta.fasta', options.tempDir + 'allFasta.fasta', options.tempDir]
	cmd = ' '.join(cmdArr)
	os.system(cmd)

	cmdArr = ['rm', options.tempDir + 'allFasta.fasta']
	cmd = ' '.join(cmdArr)
	os.system(cmd)

def getAllFeats(options):
	allFeats = {}

	f = open(options.tempDir + 'allFasta.fasta.' + str(options.kmerSize) + '.kmrs')

	for i in f:
		i = i.strip().split('\t')
		if len(i) != 2:
			continue

		if i[0] not in allFeats:
			allFeats[i[0]] = 0

	f.close()

	count = 0
	for i in sorted(allFeats):
		allFeats[i] = count
		count += 1

	return allFeats

def normalizeByTot(kmerHsh):
	tot = 0

	for i in kmerHsh:
		tot += kmerHsh[i]

	tot = float(tot)

	for i in kmerHsh:
		kmerHsh[i] /= tot

def normalizeByMarkov(kmerHsh):
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
		kmerHsh[i] /= subHsh[i]

def normalizeKMC(kmerHsh, options):
	if options == 0:
		return
	elif options == 1:
		normalizeByTot(kmerHsh)
		return
	elif options == 2:
		normalizeByMarkov(kmerHsh)


