'''
python combinedSpecModelStats.py [model dir] [spec map]
'''

from sys import argv
import glob
from sklearn.metrics import f1_score
import numpy as np
import scipy.stats

def mean_confidence_interval(data, confidence=0.95):
	a = 1.0 * np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
	return m, m-h, m+h

def getGIDSpecMap():
	f = open(argv[2])

	gidSpcHsh = {}
	for i in f:
		i = i.strip('\n').split('\t')
		gidSpcHsh[i[0]] = i[1]

	return gidSpcHsh

def parseTrue(tFile):
	f = open(tFile)

	trues = []
	for i in f:
		i = i.strip().split('\t')
		i[-1] = int(float(i[-1]))
		trues.append(i)

	f.close()

	return trues

def parsePred(pFile):
	f = open(pFile)

	preds = []
	for i in f:
		i = i.strip()
		preds.append(int(float(i)))

	f.close()

	return preds

def parseFold(tFile, gidSpcHsh, foldHsh, foldNum):
	pFile = tFile.replace('true', 'pred')

	trues = parseTrue(tFile)
	preds = parsePred(pFile)


	for i in range(0,len(trues)):
		spc = gidSpcHsh[trues[i][0]]
		ab = trues[i][1]
		if spc not in foldHsh:
			foldHsh[spc] = {'ALL': []}
		if ab not in foldHsh[spc]:
			foldHsh[spc][ab] = []
		while len(foldHsh[spc][ab]) < foldNum + 1:
			foldHsh[spc][ab].append([[],[]])
		while len(foldHsh[spc]['ALL']) < foldNum + 1:
			foldHsh[spc]['ALL'].append([[],[]])

		foldHsh[spc][ab][foldNum][0].append(trues[i][-1])
		foldHsh[spc][ab][foldNum][1].append(preds[i])
		foldHsh[spc]['ALL'][foldNum][0].append(trues[i][-1])
		foldHsh[spc]['ALL'][foldNum][1].append(preds[i])

def parseModel():
	gidSpcHsh = getGIDSpecMap()

	if argv[1][-1] != '/':
		argv[1] += '/'

	fList = glob.glob(argv[1] + '/all/*.true')

	foldHsh = {}
	for i in range(0,len(fList)):
		parseFold(fList[i], gidSpcHsh, foldHsh, i)

	return foldHsh

def getVMEME(trues, preds, prnt = False):
	vm = 0
	me = 0
	rs = 0
	ss = 0
	to = 0

	mi = 0

	for i in range(0,len(trues)):
		v = False
		m = False
		if trues[i] == 2:
			rs += 1
			if preds[i] != trues[i]:
				v = True
				vm += 1
		elif trues[i] == 0:
			ss += 1
			if preds[i] != trues[i]:
				m = True
				me += 1

		if trues[i] != preds[i]:
			mi += 1

		to += 1

	return float(vm)/rs, float(me)/ss, rs, ss, to

def getStats(foldHsh):
	statsHsh = {}
	for i in foldHsh:
		statsHsh[i] = {}
		for j in foldHsh[i]:
			# [f1, vme, me, s, r]
			statsHsh[i][j] = {
				'f1': [0,0,0],
				'vm': [0,0,0],
				'me': [0,0,0],
				'sr': [0,0,0]
			}
			f1s = []
			vms = []
			mes = []
			srs = [0,0,0]
			for k in range(0,len(foldHsh[i][j])):
				if i == 'salmonella enterica' and j == 'CIP':
					prnt = True
				f1 = f1_score(foldHsh[i][j][k][0], foldHsh[i][j][k][1], average='macro')
				vm, me, ss, rs, to = getVMEME(foldHsh[i][j][k][0], foldHsh[i][j][k][1])

				f1s.append(f1)
				vms.append(vm)
				mes.append(me)
				srs[0] += ss
				srs[1] += rs
				srs[2] += to

			statsHsh[i][j]['f1'] = list(mean_confidence_interval(f1s))
			statsHsh[i][j]['vm'] = list(mean_confidence_interval(vms))
			statsHsh[i][j]['me'] = list(mean_confidence_interval(mes))
			statsHsh[i][j]['sr'] = srs

	return statsHsh

def getABLst(statsHsh):
	abHsh = {}
	for i in statsHsh:
		for j in statsHsh[i]:
			if j not in abHsh:
				abHsh[j] = 0

	count = 0
	for i in sorted(abHsh):
		abHsh[i] = count
		count += 1

	return abHsh

def printStatsHsh(statsHsh, abHsh):
	header = ['']
	for i in sorted(abHsh):
		header += [i, '', '']
	header2 = ['Species']
	for i in range(1,len(header),3):
		header2 += ['Score', 'CI Low', 'CI High']

	keys = ['f1', 'vm', 'me', 'sr']
	for key in keys:
		header[0] = key
		print '\t'.join(header)

		if key == 'sr':
			header2 = ['Species']
			for i in range(1,len(header),3):
				header2 += ['Susceptible', 'Resistant', 'Total']

		print '\t'.join(header2)
		for spc in sorted(statsHsh):
			lArr = [''] * len(header)
			lArr[0] = spc.capitalize()
			for ab in sorted(abHsh):
				ind = abHsh[ab] * 3
				if ab in statsHsh[spc]:
					stats = statsHsh[spc][ab][key]
					lArr[ind + 1] = stats[0]
					lArr[ind + 2] = stats[1]
					lArr[ind + 3] = stats[2]
				else:
					lArr[ind + 1] = ''
					lArr[ind + 2] = ''
					lArr[ind + 3] = ''

			for i in range(0,len(lArr)):
				if isinstance(lArr[i], float):
					round(lArr[i], 3)
				lArr[i] = str(lArr[i])
			print '\t'.join(lArr)
		print ''

def main():
	foldHsh = parseModel()
	statsHsh = getStats(foldHsh)
	abHsh = getABLst(statsHsh)

	# printStatsHsh(statsHsh, abHsh)

main()