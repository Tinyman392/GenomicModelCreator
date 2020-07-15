'''
python comparePredictionsToTab.py [true tab] [pred tab] [model type]
'''

from sys import argv, stderr
from math import log
import os
from sklearn.metrics import classification_report, confusion_matrix, f1_score

def toNum(sMic):
	mic = 0
	for i in range(0,len(sMic)):
		try:
			mic = float(sMic[i:])
			break
		except:
			continue

	if '>' in sMic and '=' not in sMic:
		mic *= 2
	elif '<' in sMic and '=' not in sMic:
		mic /= 2

	mic = int(round(log(mic,2)))

	return mic

def parseTrue():
	f = open(argv[1])

	trueHsh = {}
	for i in f:
		i = i.strip().split('\t')
		i[2] = toNum(i[2])
		i[1] = i[1].split(':')[0]

		if i[0] not in trueHsh:
			trueHsh[i[0]] = {}
		if i[1] not in trueHsh[i[0]]:
			trueHsh[i[0]][i[1]] = []

		if len(trueHsh[i[0]][i[1]]) != 0:
			continue
		trueHsh[i[0]][i[1]].append(i[2])
	f.close()

	return trueHsh

def parsePred(trueHsh):
	f = open(argv[2])

	for i in f:
		i = i.strip('\n').split('\t')
		if len(i) < 3:
			continue
			
		i[0] = os.path.basename(i[0]).replace('.fasta', '')
		i[2] = int(float(i[2]))

		if i[0] not in trueHsh:
			continue
		if i[1] not in trueHsh[i[0]]:
			continue

		trueHsh[i[0]][i[1]].append(i[2])

	f.close()

def filterToAntibiotic(hsh):
	antibioHsh = {'ALL': [[], []]}
	for i in hsh:
		for j in hsh[i]:
			if len(hsh[i][j]) != 2:
				continue
			if j not in antibioHsh:
				antibioHsh[j] = [[], []]

			antibioHsh[j][0].append(hsh[i][j][0])
			antibioHsh[j][1].append(hsh[i][j][1])
			antibioHsh['ALL'][0].append(hsh[i][j][0])
			antibioHsh['ALL'][1].append(hsh[i][j][1])

	return antibioHsh

def w1(t, p, average = None):
	corr = 0
	for i in range(0,len(t)):
		if abs(t[i] - p[i]) <= 1:
			corr += 1

	return corr / len(t)

def computeScores(antibioHsh):
	scoreFunct = f1_score
	if argv[3] == 'mic' or argv[3] == 'reg':
		scoreFunct = w1

	for i in sorted(antibioHsh):
		if 1 not in antibioHsh[i][0]:
			antibioHsh[i][0].append(1)
			antibioHsh[i][1].append(1)
		score = scoreFunct(antibioHsh[i][0], antibioHsh[i][1], average = 'weighted')
		print i + '\t' + str(score)

hsh = parseTrue()
parsePred(hsh)
antibioHsh = filterToAntibiotic(hsh)
computeScores(antibioHsh)
