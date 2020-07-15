'''
python getImportantFeatures.py [model directory] [tabular file] [antibiotic] [kmc directory] 
'''

from sys import argv,stderr
from glob import glob
import xgboost as xgb
from os.path import basename

if argv[1][-1] != '/':
	argv[1] += '/'
if argv[3][-1] != '/':
	argv[3] += '/'

def err(s):
	stderr.write(s)

def parseTab():
	f = open(argv[2])

	labelHsh = {}
	for i in f:
		i = i.strip().split('\t')
		if i[0] not in labelHsh:
			labelHsh[i[0]] = {}
		
		labelHsh[i[0]][i[1]] = i[2]

	f.close()

	return labelHsh

def convToFeat(featHsh):
	f = open(argv[1] + 'model.attrOrder')

	attrHsh = {}
	for i in f:
		i = i.strip().split('\t')
		if i[1] not in featHsh:
			continue
		attrHsh[i[0]] = featHsh[i[1]]

	return attrHsh

def getImportantFeatures():
	fList = glob(argv[1] + argv[3] + '/*pkl')

	featHsh = {}
	for i in fList:
		model = xgb.Booster(model_file = i)
		hsh = model.get_score(importance_type='gain')
		for j in hsh:
			if j.replace('f', '') not in featHsh:
				featHsh[j.replace('f', '')] = 0
			featHsh[j.replace('f', '')] += hsh[j]

	return featHsh

def getKmers(fName, featHsh):
	f = open(fName)

	hsh = {}
	for i in f:
		i = i.strip().split('\t')
		if i[0] not in featHsh:
			continue
		hsh[i[0]] = 0

	f.close()

	return set(list(hsh))

def getKmerHsh(labelHsh, featHsh):
	kmerHsh = {}
	count = 0
	inc = len(labelHsh) / 50
	for i in labelHsh:
		if count >= inc:
			err('=')
			count = 0
		count += 1

		fName = glob(argv[4] + i + '.fasta.*.kmrs')
		if len(fName) == 0:
			continue
		kmerHsh[i] = getKmers(fName[0], featHsh)

	return kmerHsh

err('Parsing tab...')
labelHsh = parseTab()
err('Done\nParsing features...')
featHsh = getImportantFeatures()
err('Done\nConverting feature numbers to kmers...')
featHsh = convToFeat(featHsh)
err('Done\n')

if len(argv) < 5:
	err("Printing important feature list...")
	for i in sorted(featHsh, key = lambda x: featHsh[x], reverse = True):
		print i, featHsh[i]
	err('Done\n')
else:
	err('Reading KMC outputs...\n\t')
	kmerHsh = getKmerHsh(labelHsh, featHsh)
	err('Done\n')
