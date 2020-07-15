#! /usr/bin/env python
'''
python iterativeTest.py [fasta dir] [tabular] [temp]
'''

from sys import argv,stderr
import random
from math import log
import os
import glob
import xgboost as xgb
from sklearn.metrics import classification_report, confusion_matrix, f1_score
import numpy as np

def log2MIC(m):
	mic = 0
	for i in range(0,len(m)):
		try:
			mic = float(m[i:])
			break
		except:
			continue

	if '>' in m and '=' not in m:
		mic *= 2
	if '<' in m and '=' not in m:
		mic /= 2

	mic = int(round(log(mic, 2)))

	return mic

def parseTab():
	f = open(argv[2])

	tab = []
	gids = {}
	antibiotics = {}
	labels = {}
	for i in f:
		i = i.strip().split('\t')
		i[1] = i[1].split(':')[0]
		i[2] = log2MIC(i[2])
		tab.append(i)

		if i[0] not in gids:
			gids[i[0]] = 0

		if i[1] not in antibiotics:
			antibiotics[i[1]] = 0

		if i[2] not in labels:
			labels[i[2]] = 0

	f.close()

	count = 0
	for i in sorted(antibiotics):
		antibiotics[i] = count
		count += 1

	return tab, gids, antibiotics, labels

def runKMC(f):
	cmdArr = ['kmc.sh', '10', f, argv[3] + 'kmc/' + os.path.basename(f), argv[3]]
	cmd = ' '.join(cmdArr)

	os.system(cmd)

def parseKMC(fName):
	f = open(fName)

	kmerHsh = {}
	for i in f:
		i = i.strip().split('\t')
		i[1] = int(i[1])

		kmerHsh[i[0]] = i[1]

	f.close()

	return kmerHsh

def runKMCs(gids):
	if argv[1][-1] != '/':
		argv[1] += '/'
	if argv[3][-1] != '/':
		argv[3] += '/'

	if not os.path.exists(argv[3] + 'kmc'):
		os.mkdir(argv[3] + 'kmc/')

	for i in gids:
		if os.path.exists(argv[1] + i + '.fasta'):
			runKMC(argv[1] + i + '.fasta')

	os.system('cat ' + argv[1] + '*.fasta > ' + argv[3] + 'all.fasta')
	os.system('kmc.sh 10 ' + argv[3] + 'all.fasta ' + argv[3] + 'all.fasta ' + argv[3])

	allKmers = parseKMC(argv[3] + 'all.fasta.10.kmrs')
	count = 0
	for i in sorted(allKmers):
		allKmers[i] = count
		count += 1

	return allKmers

def buildKmerHsh(allKmers):
	fList = glob.glob(argv[3] + 'kmc/*.fasta.10.kmrs')

	kmerHsh = {}
	inc = len(fList) / 50
	count = 0
	stderr.write("Building Kmer Hsh...\n\t")
	for i in fList:
		if count >= inc:
			stderr.write('=')
			count = 0
		count += 1
		gid = os.path.basename(i).replace('.fasta.10.kmrs', '')
		genHsh = parseKMC(i)
		kmerHsh[gid] = [0] * len(allKmers)
		for j in genHsh:
			ind = allKmers[j]
			kmerHsh[gid][ind] = genHsh[j]
	stderr.write('\n')

	return kmerHsh

def shuffleTab(tab):
	tabSplit = []
	for i in range(0,10):
		tabSplit.append([])

	counts = {}
	for i in tab:
		if i[1] not in counts:
			counts[i[1]] = {}
		if i[2] not in counts[i[1]]:
			counts[i[1]][i[2]] = 0

		tabSplit[counts[i[1]][i[2]] % 10].append(i)
		counts[i[1]][i[2]] += 1

	tabShuf = []
	for i in tabSplit:
		tabShuf += i

	return tabShuf

def splitTab(tab):
	inc = len(tab) / 10
	test = tab[0:inc]
	valid = tab[inc:inc*2]
	train = tab[inc*2:]

	return train, test, valid

def toMatrix(tab, kmerHsh, antibiotics):
	X = []
	y = []
	for i in tab:
		y.append(i[2])
		antibioArr = [0]*len(antibiotics)
		antibioArr[antibiotics[i[1]]] = 1
		
		X.append(kmerHsh[i[0]] + antibioArr)

	return np.asarray(X),np.asarray(y)

def trainIter(model, params, train, kmerHsh, antibiotics, batch=1024):
	random.shuffle(train)

	count = 0
	inc = len(train) / 50
	for i in range(0,len(train),batch):
		if count >= inc:
			count = 0
			stderr.write('=')
		count += batch
		X,y = toMatrix(train[i:i+batch], kmerHsh, antibiotics)
		dMat = xgb.DMatrix(X, label = y)
		if model == None:
			model = xgb.train(params, dMat, 1, verbose_eval=False)
		else:
			model = xgb.train(params, dMat, 1, verbose_eval=False, xgb_model=model)
		del X
		del dMat
	stderr.write('\n')

	return model

def score(model, tab, kmerHsh, antibiotics, metric = f1_score):
	X,y = toMatrix(tab, kmerHsh, antibiotics)
	dMat = xgb.DMatrix(X, label = y)
	del X

	pred = model.predict(dMat)

	score = f1_score(y, pred, average = 'weighted')

	return score

tab,gids,antibiotics,labels = parseTab()
allKmers = runKMCs(gids)
kmerHsh = buildKmerHsh(allKmers)
count = 0
# for i in kmerHsh:
# 	print kmerHsh[i], i
# 	break
tab = shuffleTab(tab)

train, test, valid = splitTab(tab)

params = {'max_depth': 4, 'nthread': 20, 'num_class': len(labels), 'objective': 'multi:softmax', 'eta': 0.0625, 'silent': 1}
model = None

hist = []
for i in range(0,1000):
	stderr.write("Training epoch " + str(i) + '\n\t')
	model = trainIter(model, params, train, kmerHsh, antibiotics)

	# validation
	trScore = score(model, train, kmerHsh, antibiotics)
	teScore = score(model, test, kmerHsh, antibiotics)
	vaScore = score(model, valid, kmerHsh, antibiotics)

	print "\tEpoch training scores:", trScore, vaScore, teScore
	hist.append(vaScore)

	if len(hist) < 25:
		continue

	flag = True
	prev = ''
	for i in hist[:-25:-1]:
		if prev == '':
			prev = i
			continue

		if prev < i:
			flag = False
			break

		prev = i

	if flag:
		break


