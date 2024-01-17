import xgboost as xgb

def getSlice(dmat, num, true, options, antibiotic, tot = 10):
	sz = dmat.num_row()

	inc = float(sz) / tot

	rListTr = range(
		int(inc * ( (num + 0) % tot) ), 
		int(inc * ( (num + 1) % tot) ) 
	)
	rListVa = range(
		int( inc * ((num + 1) % tot) ), 
		int( inc * ((num + 2) % tot) )
	)
	rListTe = []

	if len(rListTr) == 0:
		rListTr = range(
			int(inc * ( (num + 0) % tot) ), 
			int( sz ) 
		)
	elif len(rListVa) == 0:
		rListVa = range(
			int( inc * ((num + 1) % tot) ), 
			int( sz )
		)

	rSetTr = set(rListTr)
	rSetVa = set(rListVa)

	for i in range(0,sz):
		if i in rSetTr or i in rSetVa:
			continue

		rListTe.append(i)

	mTr = dmat.slice(rListTr)
	mVa = dmat.slice(rListVa)
	mTe = dmat.slice(rListTe)

	f = open(options.outDir + antibiotic + '/model.' + str(num) + '.true', 'w')

	for i in rListTr:
		line = true[i]
		for j in range(0,len(line)):
			line[j] = str(line[j])
		f.write('\t'.join(line) + '\n')

	f.close()

	return mTe, mVa, mTr

def train(options, tab, labs, libsvm, antibiotic):
	params = options.parameters
	if 'max_depth' not in params:
		params['max_depth'] = options.depth
	if 'nthread' not in params:
		params['nthread'] = options.threads
	if options.row_sample != 1:
		params['subsample'] = options.row_sample
	if options.max_features != 1:
		params['colsample_bytree'] = options.max_features

	if options.classify:
		try:
			params['num_class'] = int(max(list(labs)) + 1)
		except:
			params['num_class'] = len(labs) + 1
		if 'objective' not in params:
			params['objective'] = 'multi:softmax'

	# f = open(options.outDir + '')

	# f.close()

	f = open(options.outDir + antibiotic + '/model.params', 'w')

	f.write(str(params))

	f.close()

	for i in range(0,options.foldsToRun):
		tr, va, te = getSlice(libsvm, i, tab, options, antibiotic, options.totalFolds)
		evals = [(tr,'train'),(va,'eval'),(te,'test')]

		evals_result = {}
		model = xgb.train(params, tr, options.num_rounds, evals=evals, early_stopping_rounds = options.earlyStop, evals_result = evals_result, verbose_eval = False)

		pred = model.predict(te)

		model.dump_model(options.outDir + antibiotic + '/model.' + str(i) + 'tree')
		model.save_model(options.outDir + antibiotic + '/model.' + str(i) + 'pkl')

		f = open(options.outDir + antibiotic + '/model.' + str(i) + '.pred', 'w')

		for j in pred:
			f.write(str(j) + '\n')

		f.close()

		evals = {}
		evalMet = ""
		for j in evals_result:
			count = 0
			for k in evals_result[j]:
				evalMet = k
				if count not in evals:
					evals[count] = {}
				if j not in evals[count]:
					evals[count][j] = -1

				evals[count][j] = evals_result[j][k]

		f = open(options.outDir + antibiotic + '/model.' + str(i) + '.train_hist.txt', 'w')
		for j in sorted(evals):
			f.write(str(j))
			for k in evals[j]:
				f.write('\t' + k + '\t' + str(evals[j][k]))
			f.write('\n')
		f.close()



