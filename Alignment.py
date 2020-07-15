def getAllFeats(options):
	f = open(options.alignmentFile)

	line = f.readline()

	f.close()

	line = line.split('\t')[1]

	allFeats = {}
	snpCnt = 0
	fCount = 0
	for i in range(0,len(line),5):
		for j in range(0, 5):
			feat = str(snpCnt) + '_' + str(j)
			allFeats[fCount] = feat
			fCount += 1
		snpCnt += 1

	return allFeats

def parseAlignmentFile(options):
	f = open(options.alignmentFile)

	alignmentHsh = {}
	for i in f:
		i = i.strip().split('\t')
		try:
			i[1] = list(i[1])
		except:
			i.append([])
		alignmentHsh[i[0]] = i[1]

	f.close()

	return alignmentHsh
	