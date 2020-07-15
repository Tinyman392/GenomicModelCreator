from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, ExtraTreesClassifier, ExtraTreesRegressor, BaggingClassifier, BaggingRegressor
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.svm import SVC, SVR
import numpy as np

def toMIC(micStr):
	mic = 0
	for i in range(0,len(micStr)):
		try:
			mic = float(micStr[i:])
			break
		except:
			continue

	if micStr[0] == '>' and '=' not in micStr:
		mic *= 2
	elif micStr[0] == '<' and '=' not in micStr:
		mic /= 2

	return int(round(mic,0))

def loadMatrix(options, allFeats, a):
	f = open(options.tempDir + a + '.true')

	trues = []
	for i in f:
		i = i.strip().split('\t')
		i[-1] = toMIC(i[-1])

		trues.append(i)

	f.close()

	f = open(options.tempDir + a + '.libsvm')

	y = []
	X = []
	for i in f:
		i = i.split(' ')
		if options.classify:
			i[0] = int(float(i[0]))
		else:
			i[0] = float(i[0])
		y.append(i[0])
		arr = [0] * len(allFeats)
		for j in range(1,len(i)):
			ind = int(i[j].split(':')[0])
			val = int(i[j].split(':')[1])

			arr[ind] = val

		arr = np.asarray(arr)
		X.append(arr)
	# X = np.asarray(X)

	f.close()

	return X, y, trues

def getFold(X, y, trues, weights, fold, tFolds):
	te = fold
	foldSz = len(y) / tFolds

	XTe = []
	yTe = []
	XTr = []
	yTr = []
	true = []
	foldWeights = []

	for i in range(0,tFolds):
		subX = X[i*foldSz:(i+1)*foldSz]
		suby = y[i*foldSz:(i+1)*foldSz]
		subWeights = weights[i*foldSz:(i+1)*foldSz]
		if i == te:
			XTe = subX
			yTe = suby
			true = trues[i*foldSz:(i+1)*foldSz]
		else:
			XTr = XTr + subX
			yTr = yTr + suby
			foldWeights = foldWeights + subWeights

	XTr = np.asarray(XTr)


	return np.asarray(XTr), np.asarray(yTr), np.asarray(XTe), np.asarray(yTe), true, foldWeights

def getModel(options):
	model = ''

	if options.modelType == 'RandomForest':
		if options.classify:
			model = BaggingClassifier(base_estimator=DecisionTreeClassifier(criterion = "gini", max_depth=options.depth), n_estimators=options.num_rounds, max_samples=options.row_sample, n_jobs=options.threads, max_features=options.max_features)
		else:
			model = BaggingRegressor(base_estimator=DecisionTreeRegressor(criterion = "gini", max_depth=options.depth), n_estimators=options.num_rounds, max_samples=options.row_sample, max_depth=options.depth, n_jobs=options.threads, max_features=options.max_features)
	elif options.modelType == 'ExtraTrees':
		if options.classify:
			model = ExtraTreesClassifier(n_estimators=options.num_rounds, max_samples=options.row_sample, max_depth=options.depth, n_jobs=options.threads, max_features=options.max_features)
		else:
			model = ExtraTreesRegressor(n_estimators=options.num_rounds, max_samples=options.row_sample, max_depth=options.depth, n_jobs=options.threads, max_features=options.max_features)
	elif options.modelType == 'Bagging':
		if options.classify:
			model = BaggingClassifier(base_estimator=DecisionTreeClassifier(criterion = "gini",max_depth=options.depth), n_estimators=options.num_rounds, max_samples=options.row_sample, n_jobs=options.threads)
		else:
			model = BaggingRegressor(base_estimator=DecisionTreeRegressor(criterion = "gini", max_depth=options.depth), n_estimators=options.num_rounds, max_samples=options.row_sample, max_depth=options.depth, n_jobs=options.threads)
	elif options.modelType == 'SVM':
		if options.classify:
			model = SVC()
		else:
			model = SVR()
	else:
		"Model type", options.modelType, "not supported"
		exit(1)

	return model

def train(options, tab, labs, allFeats, weights, antibiotic):
	X, y, trues = loadMatrix(options, allFeats, antibiotic)

	for i in range(0,options.foldsToRun):
		XTr, yTr, XTe, yTe, true, fWeights = getFold(X, y, trues, weights, i, options.totalFolds)
		model = getModel(options)

		f = open(options.outDir + antibiotic + '/model.' + str(i) + '.true', 'w')
		for j in true:
			for k in range(0,len(j)):
				j[k] = str(j[k])
			f.write('\t'.join(j) + '\n')
		f.close()

		model.fit(XTr, yTr, sample_weight = fWeights)
		pred = model.predict(XTe)

		f = open(options.outDir + antibiotic + '/model.' + str(i) + '.pred', 'w')
		for j in pred:
			f.write(str(j) + '\n')
		f.close()

