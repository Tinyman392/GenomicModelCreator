'''
python computeMultiSpecesCLSStats.py [model dir] [gID spec map]
'''

from sys import argv
import glob

def getGIDSpecMap():
	f = open(argv[2])

	gIDSpcHsh = {}
	for i in f:
		i = i.strip().split('\t')
		gIDSpcHsh[i[0]] = i[1]

	f.close()

	return gIDSpcHsh

def main():
	gIDSpcHsh = getGIDSpecMap()
	

main()