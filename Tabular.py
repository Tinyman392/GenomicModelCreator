from sys import stderr
from math import log
import random
import os

def convLabel(s):
	lab = 0
	for i in range(0,len(s)):
		try:
			lab = float(s[i:])
			break
		except:
			continue

	if lab == 0:
		return None

	lab = log(lab, 2)
	lab = round(lab,0)

	if '>' in s and '=' not in s:
		lab += 1
	elif '<' in s and '=' not in s:
		lab -= 1

	return lab

def loadTab(options):
	valGenHsh = {}
	if options.clustFile != '':
		f = open(options.clustFile)

		for i in f:
			i = i.strip().split('\t')[0]
			if i not in valGenHsh:
				valGenHsh[i] = 0

		f.close()

	f = open(options.tabFile)
	if options.header:
		f.readline()
		
	tab = []

	for i in f:
		i = i.strip().split('\t')

		if options.alignmentFile == '' and options.tabFeat == '' and not os.path.exists(options.fastaDir + i[0] + '.fasta'):
			stderr.write(i[0] + '\n')
			continue

		if len(i) == 1:
			continue

		if len(valGenHsh) != 0 and i[0] not in valGenHsh or '/' in i[0]:
			continue

		if len(i) == 2:
			i = [i[0], 'X:X', i[1]]
		i[1] = i[1].replace('/', '_', len(i[1])).replace(' ', '_', len(i[1]))
		if ':' not in i[1]:
			i[1] += ':'
		i[1] = i[1].split(':')
		
		if not options.enumerateClasses:
			i[-1] = i[-1].split('/')[0].split(' ')[0]
			if options.twoFold:
				i[-1] = convLabel(i[-1])
				if i[-1] is None:
					print i
					continue

		if options.noI:
			if i[-1] == '1'or i[-1] == 1:
				continue
		if options.svns:
			if i[-1] == '1':
				i[-1] = '2'
			if i[-1] == 1:
				i[-1] = 2

		tab.append([i[0], i[1][0], i[1][1], i[-1]])

	f.close()

	if options.filterClass > 0:
		counts = {}
		for i in tab:
			if i[1] not in counts:
				counts[i[1]] = {}
			if i[3] not in counts[i[1]]:
				counts[i[1]][i[3]] = 0
			counts[i[1]][i[3]] += 1

		goodAntibiotics = {}
		for i in counts:
			count = 0
			for j in counts[i]:
				if counts[i][j] > options.filterClass:
					count += 1
			if count >= 2:
				goodAntibiotics[i] = 0

		tabFilt = []
		for i in tab:
			if i[1] in goodAntibiotics:
				tabFilt.append(i)
		tab = tabFilt

	labs = {}
	gids = {}
	for i in tab:
		if i[0] not in gids:
			gids[i[0]] = 0
		if i[-1] not in labs:
			labs[i[-1]] = 0

	if options.enumerateClasses:
		count = 0
		for i in sorted(labs):
			labs[i] = count
			count += 1
	else:
		for i in sorted(labs):
			if options.classify:
				labs[i] = int(i)
			else:
				labs[i] = float(i)

	for i in range(0,len(tab)):
		tab[i][-1] = labs[tab[i][-1]]

	return tab, gids, labs

def shuffleByGid(tab):
	gidHsh = {}
	for i in tab:
		if i[0] not in gidHsh:
			gidHsh[i[0]] = i

	gids = list(gidHsh)
	shuffle(gids)
	tabShuf = []
	for i in gids:
		tabShuf.append(gidHsh[i])

	return tabShuf

def shuffleTab(tab, options):
	if options.shuffle == 'gid':
		tab = shuffleByGid(tab)
	elif options.shuffle == 'all':
		shuffle(tab)

	nFolds = options.totalFolds
	folds = []
	for i in range(0,nFolds):
		folds.append([])

	counts = {}
	for i in tab:
		if len(i) < 3:
			continue
		antibiotic = i[1]
		label = i[2]

		if antibiotic not in counts:
			counts[antibiotic] = {}
		if label not in counts[antibiotic]:
			# counts[antibiotic][label] = random.randint(0, nFolds)
			counts[antibiotic][label] = 0

		folds[counts[antibiotic][label] % nFolds].append(i)

		counts[antibiotic][label] += 1

	tab = []
	for i in folds:
		# random.shuffle(i)
		for j in i:
			tab.append(j)

	return tab

def splitByAntibiotic(tab, options):
	tabHsh = {}
	if options.individualModels:
		for i in tab:
			if i[1] not in tabHsh:
				tabHsh[i[1]] = []
			tabHsh[i[1]].append(i)
	else:
		tabHsh['all'] = tab

	return tabHsh

def getGIDs(options):
	tab, gids, labs = loadTab(options)
	gids = {}
	for i in tab:
		# print i
		if i[0] not in gids:
			gids[i[0]] = 0

	return list(gids)