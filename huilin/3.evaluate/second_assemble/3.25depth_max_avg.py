#!/bin/bash  python3
import sys
import numpy as np
# python3 '/home/huilin_hu/Documents/python笔记/3.23depth_samtools.py' '/home/huilin_hu/Assemble/reads.sort.coverage' '/home/huilin_hu/3.24avg_depth'
#此脚本求每一个contig_name的最大深度和平均深度
inputfile=sys.argv[1]
outfile=sys.argv[2]
inputFile = open(inputfile, "r")
dict1={}
dict2={}
dict={}
tatol_depth=0
for line in inputFile:
	line = line.split()
	contig_name=line[0]
	length = int(line[1])
	depth = int(line[2])
	if contig_name not in dict1:
		dict[contig_name] = []
		dict1[contig_name] = depth
	else:
		dict[contig_name].append(depth)
		dict1[contig_name] = depth + dict1[contig_name]
	#print(contig_name +"\t" + str(dict1[contig_name]))
	if contig_name not in dict2:
		dict2[contig_name] = []
	else:
		dict2[contig_name].append(length)
out=open(outfile,"w+")
avg_depth=0
for key in sorted(dict2.keys()) and sorted(dict.keys()):
	max_depth=max(dict[key])
	max_length = max(dict2[key])
	med_depth=np.median(dict[key])
	avg_depth=int(dict1[key])/int(max_length)
	#print ("all contig_name avg_depth : " + "%s\t%.2f\t%.2f\t%.2f" % (key,max_depth,avg_depth,med_depth))
	out.write("%s\t%.2f\t%.2f\t%.2f" % (key,max_depth,avg_depth,med_depth)+"\n")
out.close()
