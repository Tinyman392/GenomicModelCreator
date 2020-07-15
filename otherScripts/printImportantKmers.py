'''
python printImportantKmers.py [model directory] <plot=False>
'''

from sys import argv, stderr
from glob import glob
import xgboost as xgb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if argv[1][-1] != '/':
		argv[1] += '/'

plot = False
if len(argv) > 2:
	if argv[3][0].lower() == 't':
		plot = True

def err(s):
	stderr.write(s)

def getAttr():
	err("Getting attributes...")
	f = open(argv[1] + 'model.attrOrder')

	attrHsh = {}
	for i in f:
		i = i.strip().split('\t')
		attrHsh[i[1]] = i[0]

	f.close()
	err("Done\n")

	return attrHsh

def mergeHsh(sumHsh, hsh):
	for i in hsh:
		if i not in sumHsh:
			sumHsh[i] = 0
		sumHsh[i] += hsh[i]

def parseModel(d, sumHsh):
	fList = glob(d + '*pkl')

	err('Parsing model: ' + d + '...\n\t')
	for i in fList:
		err('=')
		mod = xgb.Booster(model_file=i)
		imp = mod.get_score(importance_type='gain')
		mergeHsh(sumHsh, imp)
	err('\n')

def isKmer(k):
	for i in k:
		i = i.lower()
		if i != 'a' and i !='c' and i != 't' and i != 'g':
			return False

	return True

def filterSumHsh(sumHsh, attrHsh):
	err("Filtering out non-kmers...\n\t")
	kmerHsh = {}
	count = 0
	inc = len(sumHsh) / 50
	for i in sumHsh:
		if count >= inc:
			err('=')
			count = 0
		count += 1
		ind = i.replace('f', '')
		kmer = attrHsh[ind]
		if not isKmer(kmer):
			continue

		kmerHsh[kmer] = sumHsh[i]
	err('\n')

	return kmerHsh

def plot(sumHsh, d):
	impArr = []
	for i in sorted(sumHsh, key = lambda x: sumHsh[x], reverse = True):
		impArr.append(sumHsh[i])

	plt.plot(range(0,len(impArr)), impArr)
	plt.title('Feature Importance Scores (Gain)')
	plt.ylabel('Feature Importance (Gain)')
	plt.xlabel('Kmer')
	plt.xlim(0,120)
	ylim = plt.ylim()
	plt.ylim(10, ylim[1])
	plt.yscale('log')

	plt.savefig(d + 'features.importance.pdf')

def printKmerHsh(kmerHsh, d):
	f = open(d + 'features.importance.tab', 'w')

	for i in sorted(kmerHsh, key = lambda x: kmerHsh[x], reverse = True):
		f.write(i + '\t' + str(kmerHsh[i]) + '\n')

	f.close()

def parseModels(attrHsh):
	dList = glob(argv[1] + '*/')
	sumHsh = {}
	for i in dList:
		parseModel(i, sumHsh)

		kmerHsh = filterSumHsh(sumHsh, attrHsh)
		plot(kmerHsh, i)
		printKmerHsh(kmerHsh, i)

	return sumHsh

attrHsh = getAttr()
sumHsh = parseModels(attrHsh)

