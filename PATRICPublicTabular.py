'''
python parsePATRICAMRTxt.py [-g genus] [-s test_standards] [-b breakpoints] [-t tabFile file]
'''

from sys import stderr
from optparse import OptionParser
from math import log
import os

def parseHeader(header):
	headHsh = {}
	header = header.strip().split('\t')
	for i in range(0,len(header)):
		headHsh[header[i]] = i

	return headHsh

def getBreakPoints(breaksFile, genus):
	f = open(breaksFile)
	f.readline()

	breakpointList = []
	for i in f:
		i = i.strip('\n').split('\t')
		g = i[0].lower()
		a = i[1].lower()
		S = float(i[2])
		R = float(i[3])

		genera = []
		for j in genus:
			if g in j:
				breakpointList.append([j, g, a, S, R])

		# i = i.strip('\n').split('\t')
		# antibiotic = i[0].lower().replace(' ', '_', len(i[0])).replace('/', '-', len(i[0]))
		# S = float(i[1])
		# R = float(i[2])
		# breakpoints[antibiotic] = [S,R]

	f.close()

	breakpoints = {}
	for i in sorted(breakpointList, key = lambda x: len(x[1]), reverse = True):
		if i[0] not in breakpoints:
			breakpoints[i[0]] = {}
		if i[2] not in breakpoints[i[0]]:
			breakpoints[i[0]][i[2]] = i[3:]

	return breakpoints

def convertMIC(mStr):
	num = 0
	for i in range(0,len(mStr)):
		try:
			num = float(mStr[i:])
			break
		except:
			continue

	num = int(log(num,2))

	if '>' in mStr and '=' not in mStr:
		num += 1
	if '<' in mStr and '=' not in mStr:
		num -= 1

	num = 2.0**num

	return num

def getSIR(mic, genus, antibiotic, breakpoints):
	if genus not in breakpoints or antibiotic not in breakpoints[genus]:
		return ''

	if mic <= breakpoints[genus][antibiotic][0]:
		return 'susceptible'
	elif mic >= breakpoints[genus][antibiotic][1]:
		return 'resistant'
	else:
		return 'intermediate'

def parseTabular(options):
	breakpoints = {}
	if options.breakpointsFile != "":
		breakpoints = getBreakPoints(options.breakpointsFile, options.genus)

	if options.tabFile == '':
		os.system("wget ftp://ftp.patricbrc.org/patric2/current_release/RELEASE_NOTES/PATRIC_genomes_AMR.txt > " + options.tempDir + '/PATRIC_genomes_AMR.txt')
		options.tabFile = options.tempDir + '/PATRIC_genomes_AMR.txt'

	f = open(options.tabFile)

	header = f.readline()
	headHsh = parseHeader(header)
	# genome_name
	# testing_standard
	# antibiotic
	# genome_id
	# resistant_phenotype
	# measurement
	
	imporperPhenotypes = {}
	improperMIC = {}

	tab = []

	for i in f:
		i = i.strip('\n').split('\t')
		genome_id = i[headHsh['genome_id']].lower()
		genome_name = i[headHsh['genome_name']].lower()
		antibiotic = i[headHsh['antibiotic']].lower()
		resistant_phenotype = i[headHsh['resistant_phenotype']].lower()
		measurement = i[headHsh['measurement']].lower().split('/')[0]
		testing_standard = i[headHsh['testing_standard']].lower()

		flag = False
		genus = ""
		for j in options.genus:
			if j in genome_name:
				flag = True
				genus = j
				break
		if len(options.genus) == 0:
			flag = True

		if not flag or (measurement == '' and testing_standard not in options.testStandard and len(options.testStandard) > 0):
			continue

		if measurement != '' and len(options.genus) > 0:
			try:
				measurement = convertMIC(measurement)
				sir = getSIR(measurement, genus, antibiotic, breakpoints)
				if sir != '':
					resistant_phenotype = sir
				elif testing_standard not in options.testStandard and len(options.testStandard) != 0:
					continue
			except:
				if measurement not in improperMIC:
					improperMIC[measurement] = 0
					stderr.write("Improper MIC: " + measurement + '\n')

		if resistant_phenotype == '':
			continue
		elif resistant_phenotype[0] == 's':
			resistant_phenotype = 0
		elif resistant_phenotype[0] == 'i':
			resistant_phenotype = 1
		elif resistant_phenotype[0] == 'r':
			resistant_phenotype = 2
		elif resistant_phenotype == 'non-susceptible':
			resistant_phenotype = 2
		elif resistant_phenotype == 'non-resistant':
			resistant_phenotype = 0
		else:
			if resistant_phenotype not in imporperPhenotypes:
				imporperPhenotypes[resistant_phenotype] = 0
				stderr.write("Improper resistant phenotype: " + resistant_phenotype + '\n')
			continue

		tab.append([genome_id, antibiotic.replace(' ', '_', len(antibiotic)).replace('/', '-', len(antibiotic)) + ':MIC', resistant_phenotype])

	f.close()

	if options.filterClass > 0:
		counts = {}
		for i in tab:
			if i[1] not in counts:
				counts[i[1]] = {}
			if i[2] not in counts[i[1]]:
				counts[i[1]][i[2]] = 0
			counts[i[1]][i[2]] += 1

		goodAntibiotics = {}
		for i in counts:
			count = 0
			for j in counts[i]:
				if len(counts[i][j]) > options.filterClass:
					count += 1
			if count >= 2:
				goodAntibiotics[i] = 0

		tabFilt = []
		for i in tab:
			if i[1] in goodAntibiotics:
				tabFilt.append(i)
		tab = tabFilt

	gids = {}
	labs = {}
	for i in tab:
		if i[0] not in gids:
			gids[i[0]] = 0
		if i[-1] not in labs:
			labs[i[-1]] = 0

	return tab, gids, labs

# parseTabular(options)
