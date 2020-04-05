#!/usr/bin/env python
#import commands
import os
import sys

program=sys.argv[0]
arguments=sys.argv[1:]
count=len(arguments)
if count !=5:
	print ("Usage: python(v3) VAF_calculate.py input_bed(file_format: chr pos-1 pos ref alt sample, sep=\"\\t\") output_features bam_dir reference_fasta n_jobs \n\nNote: 1. name of bam files are \"sample.bam\" by default. \n2.there should be a fai file under the same dir of the fasta file (samtools faidx input.fa) \n3. we did not use dbSNP AF as an feature, but you can use it to train your model if you have interest in common variants")
	sys.exit(1)
elif count==5:
	program_name = sys.argv[0]
	input_pos=sys.argv[1] #walsh.nocluster.noalt_allele_in_normal.norepeat.bed 1       1015256 1015257 A       G       Walsh
	output=sys.argv[2]
	bam_dir_tmp=sys.argv[3]
	reference_fasta=sys.argv[4]
	n_jobs=sys.argv[5]
	#sequencing_type=sys.argv[5]

from subprocess import *
from multiprocessing import Pool
import subprocess
import math
import numpy as np
import regex as re
from collections import defaultdict
from pyfaidx import Fasta
import pysam
#from scipy.stats import mannwhitneyu
base=dict()
base['A']='T'
base['T']='A'
base['G']='C'
base['C']='G'


if bam_dir_tmp.endswith('/'):
	bam_dir=bam_dir_tmp[:-1]
else:
	bam_dir=bam_dir_tmp

reference = Fasta(reference_fasta)
reference_fai =reference_fasta+".fai"
genome=reference_fai

querypos_major=defaultdict(list)
mapq_major=defaultdict(list)
baseq_major=defaultdict(list)
leftpos_major=defaultdict(list)
mismatches_major=defaultdict(list)
major_read1=dict()
major_read2=dict()
major=dict()
major_baseq=dict()
major_ids=defaultdict(list)
minor_ids=defaultdict(list)
conflict_num=dict()
baseq_major_near1b=defaultdict(list)
seqpos_major=defaultdict(list)
context1=dict()
context2=dict()
context1_count=dict()
context2_count=dict()
querypos_minor=defaultdict(list)
dp_far=defaultdict(list)
dp_near=defaultdict(list)
mapq_minor=defaultdict(list)
baseq_minor=defaultdict(list)
leftpos_minor=defaultdict(list)
mismatches_minor=defaultdict(list)
minor_read1=dict()
minor_read2=dict()
minor=dict()
minor_baseq=dict()
baseq_minor_near1b=defaultdict(list)
seqpos_minor=defaultdict(list)

#genome="/home/yd65/tools/MosaicHunter/resources/human_g1k_v37_decoy.fasta.fai"
file0=open(genome)
chr_sizes=dict()
for line in file0:
	line=line.rstrip()
	fields=line.split('\t')
	chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])
file0.close()

def process_line(line):
#file=open(input_pos)
#1       1072410 1072411 C       A       Walsh
#for line in file:
	line=line.rstrip()
	fields=line.split('\t')
	sample=fields[5]
	chr=fields[0]
	pos=int(fields[2])
	pos1=max(0,int(pos)-1)
	pos2=min(int(chr_sizes[chr]),int(pos)+1)
	major_allele=fields[3]
	minor_allele=fields[4]
	name=str(sample)+'~'+str(chr)+'~'+str(pos)+"~"+str(major_allele)+"~"+str(minor_allele)
	input_bam=bam_dir+"/"+str(sample)+".bam"
	a=pysam.AlignmentFile(input_bam, "rb")
	chrom=str(chr)
	start=int(pos)-1
	end=int(pos)
	major=dict()
	minor=dict()
	major_baseq[name]=0
	minor_baseq[name]=0
	major[name]=0
	minor[name]=0
	major_read1[name]=0
	minor_read1[name]=0
	major_read2[name]=0
	minor_read2[name]=0
	major_ids[name]=list()
	minor_ids[name]=list()
	conflict_num[name]=0
	context1_count[name]=context1_count.get(name,0)
	context2_count[name]=context2_count.get(name,0)
	mismatches_major[name]=list()
	mismatches_minor[name]=list()
	#context1[name]=reference[chrom][int(pos)-2:int(pos)+1]
	#context2[name]=(base[str(reference[chrom][int(pos)-2:int(pos)-1])]+base[str(reference[chrom][int(pos)-1:int(pos)])]+base[str(reference[chrom][int(pos):int(pos)+1])])[::-1]

	for pileupcolumn in a.pileup(chrom, start, end, max_depth=8000):
		for pileupread in pileupcolumn.pileups:
			if pileupread.indel !=0:
				continue
			try:
				querybase=pileupread.alignment.query_sequence[pileupread.query_position]
				if pileupcolumn.pos==pos-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024):
					#if sequencing_type="PE":
					#print pileupcolumn.pos
					q=pileupread.alignment.query_qualities[pileupread.query_position]
					if querybase==major_allele:
						major[name]=major.get(name,0)+1
						major_baseq[name]=major_baseq.get(name,0)+1-0.1**(float(q)/10)
						minor_baseq[name]=minor_baseq.get(name,0)+0.1**(float(q)/10)/3
						baseq_major[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
					if querybase==minor_allele:
						minor[name]=minor.get(name,0)+1
						minor_baseq[name]=minor_baseq.get(name,0)+1-0.1**(float(q)/10)
						major_baseq[name]=major_baseq.get(name,0)+0.1**(float(q)/10)/3
						baseq_minor[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
						
#      elif pileupread.query_position <1:
#       baseq_minor_near1b[name].append("end")
			except:
				continue  
#	if major_plus[name]+major_minus[name]+minor_plus[name]+minor_minus[name]>0:
#		return (name,major_plus[name],major_minus[name],minor_plus[name],minor_minus[name],major_plus[name]+major_minus[name]+minor_plus[name]+minor_minus[name], (minor_plus[name]+minor_minus[name])/(major_plus[name]+major_minus[name]+minor_plus[name]+minor_minus[name]))
#	elif major_plus[name]+major_minus[name]+minor_plus[name]+minor_minus[name]==0:
#		return (name,major_plus[name],major_minus[name],minor_plus[name],minor_minus[name],major_plus[name]+major_minus[name]+minor_plus[name]+minor_minus[name], "NA")
	if major[name]+minor[name]>0:
		return(name,chr,pos,sample,major_allele, minor_allele, major[name],minor[name],major_baseq[name],minor_baseq[name],','.join(str(x) for x in baseq_major[name])+",",','.join(str(x) for x in baseq_minor[name])+",",major[name]+minor[name],minor[name]/(major[name]+minor[name]))
	elif major[name]+minor[name]==0:
		return(name,chr,pos,sample,major_allele, minor_allele, major[name],minor[name],major_baseq[name],minor_baseq[name],','.join(str(x) for x in baseq_major[name])+",",','.join(str(x) for x in baseq_minor[name])+",",0,"NA")
			
			
fo=open(output,"w")
#header='id major_plus major_minor minor_plus minor_minus depth AF '.split()
header='id sample major_allele minor_allele major_count minor_count major_count_baseq minor_count_baseq major_baseQs minor_baseQs depth AF'.split()
print ('\t'.join(header),file=fo)

if __name__ == "__main__":
	pool = Pool(processes=int(n_jobs))
	with open(input_pos) as source_file:
	# chunk the work into batches of 4 lines at a time
		#pool.map(process_line, source_file,1)
		result=pool.map(process_line, source_file,1)
		for atuple in result:
			try:
				print ('\t'.join(str(x) for x in atuple),file=fo)
			except:
				continue


fo.close()




