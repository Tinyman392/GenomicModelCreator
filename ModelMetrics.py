import numpy as np
import scipy.stats
from glob import glob
from math import log
from sklearn.metrics import classification_report, confusion_matrix, f1_score, r2_score, recall_score, precision_score

def mean_confidence_interval(data, confidence=0.95):
	a = 1.0 * np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
	return m, m-h, m+h

def parseTrue(fName, toFloat = False, roundPred = False):
	f = open(fName)

	true = []
	for i in f:
		i = i.strip('\n').split('\t')
		antibiotic = i[1].split(':')[0]
		lab = i[-1]

		if toFloat:
			lab = float(i[-1])
		if roundPred:
			lab = int(round(float(i[-1]), 0))

		true.append([antibiotic, lab])

	f.close()

	return true

def parsePred(fName, toFloat = False, roundPred = False):
	f = open(fName)

	pred = []
	for i in f:
		i = i.strip()

		if toFloat:
			i = float(i)
			if roundPred:
				i = int(round(i, 0))

		pred.append(i)

	f.close()

	return pred

def computeR2(options, a):
	predFs = sorted(glob(options.outDir + a + '/*.pred'))
	trueFs = sorted(glob(options.outDir + a + '/*.true'))

	r2s = []
	for i in range(0,len(predFs)):
		true = parseTrue(trueFs[i], True)
		pred = parsePred(predFs[i], True, False)

		for j in range(0,len(true)):
			true[j] = true[j][-1]

		r2 = r2_score(true, pred)
		r2s.append(r2)

	conf = list(mean_confidence_interval(r2s))
	for i in range(0,len(conf)):
		conf[i] = str(conf[i])

	fName = options.outDir + a + '/r2.tab'
	f = open(fName, 'w')

	f.write('\t'.join(conf) + '\n')

	f.close()

def computeAcc(options, func, fName, a):
	predFs = sorted(glob(options.outDir + a + '/*.pred'))
	trueFs = sorted(glob(options.outDir + a + '/*.true'))

	accs = []
	antibiotics = {'all':0}

	for i in range(0,len(predFs)):
		true = parseTrue(trueFs[i], True)
		pred = parsePred(predFs[i], True, True)

		accHsh = {'all': [0, 0]}
		for j in range(0,len(true)):
			t = true[j][1]
			antibiotic = true[j][0]
			p = pred[j]

			if antibiotic not in antibiotics:
				antibiotics[antibiotic] = 0

			if antibiotic not in accHsh:
				accHsh[antibiotic] = [0,0]
			if func(t, p):
				accHsh[antibiotic][0] += 1
				accHsh['all'][0] += 1
			accHsh[antibiotic][1] += 1
			accHsh['all'][1] += 1

		for j in accHsh:
			accHsh[j][0] = float(accHsh[j][0]) / accHsh[j][1]

		accs.append(accHsh)

	f = open(fName, 'w')

	f.write('Antibiotic\tAcc\tConf_low\tConf_high\tCount\n')
	for i in sorted(antibiotics):
		antibioticAcc = []
		antibioticCnt = 0

		for j in accs:
			if i in j:
				antibioticAcc.append(j[i][0])
				antibioticCnt += j[i][1]

		conf = list(mean_confidence_interval(antibioticAcc))

		arr = [i] + conf + [antibioticCnt]
		for j in range(0,len(arr)):
			arr[j] = str(arr[j])
		f.write('\t'.join(arr) + '\n')

	f.close()

def equals(a, b):
	return a == b

def w1(a, b):
	return abs(a - b) <= 1

def getRawAcc(options, a):
	computeAcc(options, equals, options.outDir + a + '/raw_acc.tab', a)

def getW1Acc(options, a):
	computeAcc(options, w1, options.outDir + a + '/w1_acc.tab', a)

def getClassificationReport(options, labs, a):
	pass

	# rLabs = {}
	# for i in labs:
	# 	rLabs[labs[i]] = i

	# predFs = sorted(glob(options.outDir + a + '/*.pred'))
	# trueFs = sorted(glob(options.outDir + a + '/*.true'))

	# clVals = {'all':[[],[]]}

	# for i in range(0,len(predFs)):
	# 	true = parseTrue(trueFs[i], True)
	# 	pred = parsePred(predFs[i], True, True)

	# 	for j in range(0,len(true)):
	# 		t = true[j][-1]
	# 		p = pred[j]

	# 		tLab = rlabs[t]

	# 		if tLab not in clVals:
	# 			clVals[tLab] = [[],[]]

	# 		clVals[tLab][0].append(t)
	# 		clVals[tLab][1].append(p)

	# 		clVals['all'][0].append(t)
	# 		clVals['all'][1].append(p)

	# clStats = {}
	# stats = [f1_score, precision_score, recall_score]
	# for i in clVals:
	# 	if i not in clStats:
	# 		clStats[i] = [[], [], []]
	# 	f = f1_score(clStats[i][0], clStats[i][1], average='macro')
	# 	p = precision_score(clStats[i][0], clStats[i][1], average='macro')
	# 	r = recall_score(clStats[i][0], clStats[i][1], average='macro')

	# 	clStats[i][0].append(p)
	# 	clStats[i][1].append(r)
	# 	clStats[i][2].append(f)

	# for i in clStats:
	# 	for j in range(0,len(clStats[i])):
	# 		clStats[i][j] = list(mean_confidence_interval(clStats[i][j]))

def getAntibiotics(options, a):
	antibiotics = {}
	trueFs = sorted(glob(options.outDir + a + '/*.true'))

	for i in trueFs:
		f = open(i)

		for j in f:
			j = j.strip().split('\t')[1].split(':')[0]
			if j not in antibiotics:
				antibiotics[j] = 0

		f.close()

	return antibiotics

def getBreaks(options, a, labs):
	breaks = {}
	if options.breakpointsFile == '':
		antibiotics = getAntibiotics(options, a)
		for i in antibiotics:
			breaks[i] = [1, 4]
	else:
		f = open(options.breakpointsFile)
		for i in f:
			if i[0] == '#':
				continue
			i = i.strip('\n').split('\t')
			try:
				if len(i) == 4:
					i[0] = i[0].replace('/','_',len(i[0])).replace(' ','_',len(i[0]))
					breaks[i[0]] = [float(i[1]), float(i[2])]
				if len(i) > 4:
					i[1] = i[1].replace('/','_',len(i[1])).replace(' ','_',len(i[1]))
					if options.genus in i[0]:
						breaks[i[1]] = [float(i[2]), float(i[3])]
			except:
				continue

	for i in sorted(breaks):
		for j in range(0,len(breaks[i])):
			breaks[i][j] = int(round(log(breaks[i][j], 2)))
		# print i, breaks[i]

	return breaks

def toSIR(m, breaks, antibiotic):
	sir = 'I'
	if antibiotic not in breaks:
		return ''

	if m >= breaks[antibiotic][1]:
		return 'R'
	if m <= breaks[antibiotic][0]:
		return 'S'

	return sir

def getVMEME(options, a):
	predFs = sorted(glob(options.outDir + a + '/*.pred'))
	trueFs = sorted(glob(options.outDir + a + '/*.true'))

	true = parseTrue(trueFs[0], True)
	labs = {}
	for i in true:
		if i[-1] not in labs:
			labs[i[-1]] = 0
	labs = list(labs)

	breaks = getBreaks(options, a, labs)

	accs = []
	antibiotics = {'all':0}

	for i in range(0,len(predFs)):
		true = parseTrue(trueFs[i], True, True)
		pred = parsePred(predFs[i], True, True)

		for j in range(0,len(true)):
			true[j][-1] = int(true[j][-1])
			pred[j] = int(pred[j])

		# accHsh[antibiotic] -> [VME, R, ME, S]
		accHsh = {'all': [0, 0, 0, 0]}
		for j in range(0,len(true)):
			t = true[j][1]
			antibiotic = true[j][0]
			p = pred[j]

			tSIR = toSIR(t, breaks, antibiotic)
			pSIR = toSIR(p, breaks, antibiotic)

			if antibiotic not in antibiotics:
				antibiotics[antibiotic] = 0

			if antibiotic not in accHsh:
				accHsh[antibiotic] = [0,0,0,0]

			if tSIR == 'R':
				if pSIR == 'S':
					accHsh['all'][0] += 1
					accHsh[antibiotic][0] += 1
				accHsh['all'][1] += 1
				accHsh[antibiotic][1] += 1
			if tSIR == 'S':
				if pSIR == 'R':
					accHsh['all'][2] += 1
					accHsh[antibiotic][2] += 1
				accHsh['all'][3] += 1
				accHsh[antibiotic][3] += 1

		for j in accHsh:
			if accHsh[j][1] == 0:
				accHsh[j][0] = -1
			else:
				accHsh[j][0] = float(accHsh[j][0]) / accHsh[j][1]
			if accHsh[j][3] == 0:
				accHsh[j][2] = -1
			else:
				accHsh[j][2] = float(accHsh[j][2]) / accHsh[j][3]

		accs.append(accHsh)

	f = open(options.outDir + a + '/VMEME.tab', 'w')
	f.write('Antibiotic\tVME_Acc\tVME_Conf_low\tVME_Conf_high\tME_Acc\tME_Conf_low\tME_Conf_high\tRes_Count\tSus_Count\n')

	for i in sorted(antibiotics):
		vmeAcc = []
		vmeCnt = 0
		meAcc = []
		meCnt = 0

		for j in accs:
			if i in j:
				if j[i][0] != -1:
					vmeAcc.append(j[i][0])
				vmeCnt += j[i][1]
				if j[i][2] != -1:
					meAcc.append(j[i][2])
				meCnt += j[i][3]

		vmeConf = list(mean_confidence_interval(vmeAcc))
		meConf = list(mean_confidence_interval(meAcc))

		arr = [i] + vmeConf + meConf + [vmeCnt, meCnt]
		for j in range(0,len(arr)):
			arr[j] = str(arr[j])
		f.write('\t'.join(arr) + '\n')

	f.close()

def catPredTrue(options, a):
	predFs = sorted(glob(options.outDir + a + '/*.pred'))
	trueFs = sorted(glob(options.outDir + a + '/*.true'))

	trues = []
	preds = []
	for i in range(0,len(predFs)):
		true = parseTrue(trueFs[i], True)
		pred = parsePred(predFs[i], True, True)

		trues = trues + true
		preds = preds + pred

	return trues, preds

def getLabArr(labs):
	labArr = []
	try:
		labArr = [''] * (int(max(list(labs))) + 1)
	except:
		print labs
		labArr = [''] * (len(labs) + 1)
	for i in labs:
		try:
			labArr[labs[i]] = i
		except:
			labArr[int(labs[i])] = int(i)

	return labArr

def classificationReport(options, labs, a):
	trues, preds = catPredTrue(options, a)
	labArr = getLabArr(labs)

	dLst = []
	for i in range(0,len(labArr)):
		if labArr[i] == '':
			dLst.append(i)
	dLst = dLst[::-1]

	for i in dLst:
		del labArr[i]

	trues = convArrToLabs(trues, labArr)
	preds = convArrToLabs(preds, labArr)

	report = classification_report(trues, preds, labels = labArr)

	f = open(options.outDir + a + '/classification_report.txt', 'w')
	f.write(report)
	f.close()

def convArrToLabs(arr, labArr):
	conv = []
	for i in arr:
		if type(i) == list:
			i = i[-1]
		
		i = int(i)

		conv.append(labArr[i])

	return conv

def confusionMatrix(options, labs, a):
	trues, preds = catPredTrue(options, a)
	labArr = getLabArr(labs)

	trues = convArrToLabs(trues, labArr)
	preds = convArrToLabs(preds, labArr)

	for i in labArr:
		if i not in trues:
			trues.append(i)
			preds.append(i)

	confMat = confusion_matrix(trues, preds, labArr)
	for i in range(0,len(confMat)):
		for j in range(0,len(confMat[i])):
			confMat[i][j] = str(confMat[i][j])

	f = open(options.outDir + a + '/confusion_matrix.tab', 'w')

	for i in labArr:
		f.write('\t' + str(i))
	f.write('\n')

	for i in range(0,len(confMat)):
		f.write(str(labArr[i]))
		for j in confMat[i]:
			f.write('\t' + str(j))
		f.write('\n')

	f.close()

def f1(options, labs, a):
	predFs = sorted(glob(options.outDir + a + '/*.pred'))
	trueFs = sorted(glob(options.outDir + a + '/*.true'))

	fHsh = {}

	# print labs

	for i in range(0,len(predFs)):
		pred = parsePred(predFs[i], True, True)
		true = parseTrue(trueFs[i], True, True)

		aHsh = {'ALL': [[], []]}
		for j in range(0,len(true)):
			antibiotic = true[j][0]
			t = true[j][1]
			p = pred[j]

			if antibiotic not in aHsh:
				aHsh[antibiotic] = [[],[]]

			aHsh[antibiotic][0].append(t)
			aHsh[antibiotic][1].append(p)
			aHsh['ALL'][0].append(t)
			aHsh['ALL'][1].append(p)

		for j in aHsh:
			# delLst = []
			# for k in range(0,len(aHsh[j][0])):
			# 	if aHsh[j][1][k] == 1:# or aHsh[j][0][k] == 1:
			# 		delLst.append(k)
			# delLst = delLst[::-1]
			# for k in delLst:
			# 	del aHsh[j][0][k]
			# 	del aHsh[j][1][k]

			f = f1_score(aHsh[j][0], aHsh[j][1], average='macro')
			if j not in fHsh:
				fHsh[j] = []
			fHsh[j].append(f)

		# break

	f = open(options.outDir + a + '/f1.tab', 'w')

	for i in sorted(fHsh):
		conf = mean_confidence_interval(fHsh[i])
		conf = list(conf)
		arr = [i] + conf

		for j in range(0,len(arr)):
			arr[j] = str(arr[j])

		f.write('\t'.join(arr) + '\n')

	f.close()

		

