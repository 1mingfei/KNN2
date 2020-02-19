import numpy as np

def clusterRatio(fname, topN):
	'''
	fname : XYZ format file name
	topN : get Mg/Zn ratio for top N large clusters
	'''
	with open(fname, 'r') as fin:
		lines= fin.readlines()
	tot_num=int(lines[0].split()[0])
	NMg, NZn = 0, 0
	for i in range(2, 2 + tot_num):
		if (int(lines[i].split()[4]) < (topN + 1)):
			if lines[i].split()[3] == "Mg" :
				NMg += 1
			else :
				NZn += 1
	return(float(NMg / NZn))

def clusterRatioSize(fname, topN):
	'''
	fname : XYZ format file name
	topN : get Mg/Zn ratio for top N large clusters
	'''
	with open(fname, 'r') as fin:
		lines= fin.readlines()
	tot_num=int(lines[0].split()[0])
	res = []
	for i in range(1, topN + 1):
		NMg, NZn, NSize = 0, 0, 0
		for j in range(2, 2 + tot_num):
			if (int(lines[j].split()[4]) == i):
				NSize += 1
				if lines[j].split()[3] == "Mg" :
					NMg += 1
				else :
					NZn += 1
		res.append([NSize, float(NMg / NZn)])
		updateRatioSizeMap([NSize, float(NMg / NZn)])

	return res

ratioSizeMap = {}
def updateRatioSizeMap(lst):
	if (str(lst[0]) not in ratioSizeMap):
		ratioSizeMap[str(lst[0])] = [lst[1], 1]
	else :
		tmpVal = float(ratioSizeMap[str(lst[0])][0] * ratioSizeMap[str(lst[0])][1])
		tmpVal += float(lst[1])
		tmpVal /= (ratioSizeMap[str(lst[0])][1] + 1)
		ratioSizeMap[str(lst[0])] = [tmpVal, (ratioSizeMap[str(lst[0])][1]+1)]
	return

for i in range(1, 382):
	res = clusterRatioSize(str(i) + ".xyz", 6)
	#print(i, res)


ratioSizeMapSimple = {}
for key in ratioSizeMap:
	print(key, '->', ratioSizeMap[key])
	ratioSizeMapSimple[key] = ratioSizeMap[key][0]

for key in ratioSizeMapSimple:
	print(key, '->', ratioSizeMapSimple[key])

lists = sorted(ratioSizeMapSimple.items())

import matplotlib.pyplot as plt
x, y = zip(*lists)
plt.plot(x, y)
plt.xlabel("cluster size")
plt.ylabel("Mg/Zn ratio")
plt.show()