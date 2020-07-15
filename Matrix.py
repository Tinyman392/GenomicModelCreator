import xgboost as xgb

def tabToMat(tab, kmerHsh, allFeats):
	X = []
	y = []
	for i in tab:
		arr = [0]*len(allFeats)
		gid = i[0]
		arr[allFeats[i[1]]] = 1
		for j in kmerHsh[gid]:
			arr[allFeats[j]] = kmerHsh[gid][j]
		X.append(arr)
		y.append(i[2])

	return X,y

def makeMatrix(options, tab, kmerHsh, allFeats, fold=0):
	antibioticHsh = getAntibioticHsh(tab)

	trT = []
	teT = []
	vaT = []

	te = fold % options.totalFolds
	va = (fold + 1) % options.totalFolds
	foldSz = len(i) / options.totalFolds
	for i in range(0, options.totalFolds):
		fold = tab[i:i+foldSz]

		if i == te:
			teT = teT + fold
		elif i == va:
			vaT = vaT + fold
		else:
			trT = trT + fold

	trX, trY = tabToMat(trT, kmerHsh, allFeats)
	teX, teY = tabToMat(teT, kmerHsh, allFeats)
	vaX, vay = tabToMat(vaT, kmerHsh, allFeats)

	tr = xgb.DMatrix(trX, label=trY)
	te = xgb.DMatrix(teX, label=teY)
	va = xgb.DMatrix(vaX, label=vaY)

	return tr, te, va

