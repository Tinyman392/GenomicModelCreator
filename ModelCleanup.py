import os
import shutil

def moveFiles(options, gids, allFeats, antibiotic, labs, tabHsh):
	if options.cleanup:
		for i in tabHsh:
			# shutil.rmtree(options.tempDir + 'kmc/')
			if os.path.exists(options.kmcDir) and '/dev/shm/' in options.kmcDir:
				shutil.rmtree(options.kmcDir)
			if os.path.exists(options.tempDir + i + '.libsvm'):
				os.remove(options.tempDir + i + '.libsvm')


	f = open(options.outDir + 'model.params', 'w')
	f.write(str(options))
	f.close()

	f = open(options.outDir + 'model.genomes.list', 'w')
	for i in sorted(gids):
		f.write(i + '\n')
	f.close()

	f = open(options.outDir + 'model.attrOrder', 'w')
	for i in sorted(allFeats):
		f.write(str(i) + '\t' + str(allFeats[i]) + '\n')
	f.close()

	f = open(options.outDir + 'model.labels.map', 'w')
	for i in sorted(labs, key = lambda x: labs[x]):
		f.write(str(i) + '\t' + str(labs[i]) + '\n')
	f.close()

	
