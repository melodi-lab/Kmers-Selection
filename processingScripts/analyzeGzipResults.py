import os
import gzip
import sys
#from matplotlib.pyplot import *

if len(sys.argv) != 2:
	print "Require one argument. Exit"
	exit(-1)

fileName = sys.argv[1]


#fileName = 'hg38.width=10000.offset=10000.k=3.txt'
#fileName = 'hg38.width=1000000.offset=1000000.k=3.txt'
refOrder = fileName + ".ranked.ref"
#refOrder = "${fileName}.ranked.ref"
#length = `less $refOrder | wc -l`
refLines = open(refOrder,'r').readlines()
#totLines = len(fileLines)
textLines = open(fileName,'r').readlines()

# read the references in order

refList = {};
index = 0;
for refline in refLines:
	#print refline
	refline = refline.strip('\n')
	#print refline
	refList[refline] = index
	index = index + 1
	#print refList[parts[0]]
	#print parts[0]


def getKey(line):
	parts = line.split()
	return int(refList[parts[0]])


txtRefList = []
txtOrder = []
textLines = textLines[1:]
textLines = sorted(textLines, key=getKey)

NBatches = 10
compressionRatioList = []

for numBatches in range(NBatches):
	stopLines = (numBatches+1)*len(textLines)/NBatches
	newFile = fileName + '.B' + str(numBatches)
	f = open(newFile, 'w')
	for index in range(len(textLines)):
		f.write(textLines[index])		
		index = index + 1
		if index >= stopLines:
			break;
	f.close()
	statinfo = os.stat(newFile)
	origSize = statinfo.st_size
	os.system('gzip ' + newFile);
	"""
	with gzip.open(newFile+'.gz', 'wb') as f:
		#for index in range(stopLines):
		#f.write(textLines[index])
		f.write("\n".join(textLines[0:stopLines]))
		#f.write(textLines[0:stopLines].join("\n"))
	"""
	compressedSize = os.stat(newFile+'.gz').st_size	
	print "Original file size is " + str(origSize) +  "; Compressed file size is " + str(compressedSize) + "; The compression ratio is " + str(float(origSize)/compressedSize)
	compressionRatioList.append(float(origSize)/compressedSize)


f = open(fileName+'.compressedResult.txt', 'wt')
for i in compressionRatioList:
	f.write(str(i) + '\n')
f.close()
#plot(range(len(compressionRatioList)+1), compressionRatioList)
#show()

"""		
for index in range(len(textLine)):
		
	index = index + 1


#textLines = [ textLines[i]  for i in txtOrder ]
#textLines = textLines[txtOrder]

debugFile = fileName + '.debug'
f = open(debugFile, 'w')
for line in textLines:
	f.write(line)
	
f.close()

#print totLines

"""
