'''
python combinedSpecModelMICStats.py [model dir] [spc map] <make mat = False>
'''

from sys import argv
import glob
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

def getW1(arrs):
	trues = arrs[0]
	preds = arrs[1]

	sArr = []
	for i in range(0,len(trues)):
		t = trues[i]
		p = preds[i]

		if abs(t - p) <= 1:
			sArr.append(1)
		else:
			sArr.append(0)

	return float(sum(sArr)) / len(sArr), len(sArr)

def computeW1s(foldHsh):
	for i in foldHsh:
		for j in foldHsh[i]:
			for k in range(0,len(foldHsh[i][j])):
				foldHsh[i][j][k] = list(getW1(foldHsh[i][j][k]))

	for i in foldHsh:
		for j in foldHsh[i]:
			arr = [[],[]]
			for k in range(0,len(foldHsh[i][j])):
				arr[0].append(foldHsh[i][j][k][0])
				arr[1].append(foldHsh[i][j][k][1])
			arr[1] = sum(arr[1])
			arr[0] = list(mean_confidence_interval(arr[0]))
			arr = arr[0] + [arr[1]]
			foldHsh[i][j] = arr

def main():
	foldHsh = parseModel()
	computeW1s(foldHsh)

	makeMat = True
	if len(argv) > 3 and argv[3][0].lower() == 'f':
		makeMat = False

	if not makeMat:
		for i in sorted(foldHsh):
			for j in sorted(foldHsh[i]):
				arr = [i,j] + foldHsh[i][j]
				for k in range(0,len(arr)):
					arr[k] = str(arr[k])
				print '\t'.join(arr)
	else:
		ab = {}
		for i in sorted(foldHsh):
			for j in sorted(foldHsh[i]):
				if j not in ab:
					ab[j] = 0

		count = 0
		for i in sorted(ab):
			ab[i] = count
			count += 1

		arr = ['species']
		for i in sorted(ab):
			arr.append(i)
		print '\t'.join(arr)
		for i in sorted(foldHsh):
			arr = ['']*len(ab)
			for j in sorted(foldHsh[i]):
				val = str(foldHsh[i][j][0])
				ind = ab[j]
				arr[ind] = val
			arr = [i] + arr
			print '\t'.join(arr)


if __name__ == '__main__':
	main()