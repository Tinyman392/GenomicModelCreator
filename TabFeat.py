def getAllFeats(options):
	f = open(options.tabFeat)

	allFeats = {}
	ln = f.readline().split('\t')
	for i in range(0,len(ln)-1):
		allFeats[i] = i

	f.close()

	return allFeats

def parseFeats(options):
	f = open(options.tabFeat)

	tabFeatHsh = {}
	for i in f:
		i = i.strip('\n').split('\t')
		gid = i[0]
		feats = i[1:]
		for j in range(0,len(feats)):
			feats[j] = float(feats[j])
		tabFeatHsh[gid] = feats

	f.close()

	return tabFeatHsh