#! /usr/bin/env python
from __future__ import print_function
import sys
import math
import numpy as np
import optparse
import re
import random

def shuffled(x):
    y = x[:]
    random.shuffle(y)
    return y

# Function to find position of position information in FORMAT field 
def findPosInFORMATStr(infoStr, str2ID):
	RETURNPOS = 0
	info = infoStr.split(":")
	for i in info:
		if i == str2ID:
			break
		RETURNPOS = RETURNPOS+1
	return(RETURNPOS)
	

parser = optparse.OptionParser()

parser.add_option('-q', action="store", dest="QUALThreshold", type="float") # NO LONGER USED !
parser.add_option('-m', '--min_cov',action="store", dest="minCoverageThreshold", type="int")
parser.add_option('-M', '--max_cov',action="store", dest="maxCoverageThreshold", type="int")
parser.add_option('-f', '--vcf_file',action="store", dest="namefp", type="string")
parser.add_option('-c', '--CI99pc_cov',action="store",dest="namefile99CI",type="string")

(options, args) = parser.parse_args()

minCoverageThresholdIni = int(options.minCoverageThreshold)
maxCoverageThresholdIni = int(options.maxCoverageThreshold)

##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
GQThreshold = float(options.QUALThreshold)

namefile = str(options.namefp)
CIfile=str(options.namefile99CI)
fp = open(namefile)
CI = open(CIfile)

############## MAIN ############
nameSeq = []
nameInd = []
sequenceDNA = ""


seqName = ""
countSNP = 0
countSite = 0
NumberOfInd = ""
sites = []

pos = 0
oldpos = 0 

# read 99% confidence intervals COVERAGE

dicoCOV={} # store quantiles in this dictionary
for line in CI:
    line = line.strip("\n")
    line = line.replace('\"','')
    #print(line)
    if "V1" not in line: 
        splittedline = line.split()
        ID = str(splittedline[0])
        dicoCOV[ID]=splittedline
        newline = str(ID) + "\t" + line
        print(newline)

print(dicoCOV)
CI.close()
# read the vcf file

indexpos=0
for line in fp:
	line = line.rstrip()
	
	# count number of individuals and store names in array
	if line[0:6] == "#CHROM" and len(nameInd) == 0:
		
		arrline = line.split()
		
		for i in arrline[9:]:
		#	indseq1 = i+".1"
		#	indseq2 = i+".2"
		#	nameSeq.append(indseq1)
		#	nameSeq.append(indseq2)
			nameInd.append(i)
		#NumberOfInd = len(nameSeq)
		print(nameInd)
		print(nameInd[0])
	if line[0] != "#":
		indexpos+=1
		if indexpos == 1: # first loop initiate sequences
			ploidy={} # store ploidy in this dictionary
			numberindividual=0
			firstline = line.split()
			for genotype in firstline[9:]:
				genotypeGT=genotype.split(":")[0]
				if genotypeGT.count('/') == 1:	
					print("ploidy=2 ",nameInd[numberindividual])
					ploidy[nameInd[numberindividual]]=2
					print(ploidy[nameInd[numberindividual]])
					indseq1 = nameInd[numberindividual]+".1"
					indseq2 = nameInd[numberindividual]+".2"
					nameSeq.append(indseq1)
					nameSeq.append(indseq2)
				elif genotypeGT.count('/') == 2:
					print("ploidy=3 ",nameInd[numberindividual])
					ploidy[nameInd[numberindividual]]=3
					print(ploidy[nameInd[numberindividual]])
					indseq1 = nameInd[numberindividual]+".1"
					indseq2 = nameInd[numberindividual]+".2"
					indseq3 = nameInd[numberindividual]+".3"
					nameSeq.append(indseq1)
					nameSeq.append(indseq2)
					nameSeq.append(indseq3)
				elif genotypeGT.count('/') == 3:
					print("ploidy=4 ",nameInd[numberindividual])
					ploidy[nameInd[numberindividual]]=4
					print(ploidy[nameInd[numberindividual]])
					indseq1 = nameInd[numberindividual]+".1"
					indseq2 = nameInd[numberindividual]+".2"
					indseq3 = nameInd[numberindividual]+".3"
					indseq4 = nameInd[numberindividual]+".4"
					nameSeq.append(indseq1)
					nameSeq.append(indseq2)
					nameSeq.append(indseq3)
					nameSeq.append(indseq4)
				numberindividual+=1
			NumberOfInd = len(nameSeq)
			print(NumberOfInd)
		arrline = line.split()
		chromosome = arrline[0]
		pos = int(arrline[1])
		REF = str(arrline[3])
		ALT = str(arrline[4])
		site = []
		
		if (pos % 50000) == 0:
			print("Positions "+str(pos)+" of "+chromosome)

		if chromosome != seqName and seqName != "": # print seq when new chromosome found
			print("Print  " + seqName)
			fasta = open(seqName+".fst", "w")
			sites = np.ravel(sites)
			print(sites)
			sequenceDNA = np.reshape(sites, newshape=(len(sites)/NumberOfInd, NumberOfInd))
			for i, seqName in enumerate(nameSeq):
					print(">"+seqName, file=fasta)
					print(''.join(map(str, sequenceDNA[1:,i])), file=fasta)
					
			seqName = chromosome
			sites = []
			
		if seqName == "":
			seqName = chromosome
		
		
		### MODIF PROPOSER PAR EMERIC START###############################
		## This was at the end before
		if oldpos+1 < pos: # if current site is not exactly one position after oldpos site put "N"
			numberSite = int(pos-oldpos)
			for j in range(1,numberSite):
				site = []
				#print(str(pos))
				for i in range(0, len(nameSeq)): # already take care of the ploidy since it is for all nameSeq (not only nameInd)
					site.append("N")
				#if len(site) != 40:
				#	print(str(pos))
				sites = sites+site
		### MODIF PROPOSER PAR EMERIC STOP###############################
		
		#print(len(sites)," ",len(REF))
		if len(REF) > 1 or len(ALT) > 1 or ALT == "*" :  # exlcude indels
			#print("site:",len(sites)/NumberOfInd," nameseq=",len(nameSeq))
			for i in range(0, len(nameSeq)): # already take care of the ploidy since it is for all nameSeq (not only nameInd)
				site.append("N")
			#continue
			#print("excluded ",len(sites)," ",len(REF))
				
		else:

			covPos = findPosInFORMATStr(arrline[8], "DP") # find coverage position in INFO field

			if ALT == ".": # if no ALT allele ( = monorphic site)	
				if arrline[6] == "LowQual" or arrline[6] == "PASS" or arrline[6] == "FAILED_MQ;LowQual": # exclude LowQual for HaplotypeCaller and GenotypeGVCFs SNP - I add here the case of weird PASS sites with ALT == "." in GATK4
					for i in range(0, len(nameSeq)):
						site.append("N")
					#continue
				else:
					numberindividual=0
					for i in range(0, len(nameInd)):
						ind = arrline[i+9]
						arrInd = ind.split(":")
						minCoverageThreshold = minCoverageThresholdIni
                                        	maxCoverageThreshold = maxCoverageThresholdIni
                                        	#print(str(len(arrInd)) + " " + str(covPos))
						if len(arrInd) < int(covPos+1): # if individual doesn't have DP field
							if ploidy[nameInd[numberindividual]]==2:
							#print("here")
								site.append("N")
								site.append("N")
							elif ploidy[nameInd[numberindividual]]==3:
								site.append("N")
								site.append("N")
								site.append("N")
							elif ploidy[nameInd[numberindividual]]==4:
								site.append("N")
								site.append("N")
								site.append("N")
								site.append("N")
							else:
								break
						else:
							cov = arrInd[covPos]
							if cov == ".":
								if ploidy[nameInd[numberindividual]]==2:
								#print("here")
									site.append("N")
									site.append("N")	
									continue
								elif ploidy[nameInd[numberindividual]]==3:
									site.append("N")
									site.append("N")
									site.append("N")
									continue
								elif ploidy[nameInd[numberindividual]]==4:
									site.append("N")
									site.append("N")
									site.append("N")
									site.append("N")
									continue
								else:
									break
							cov = int(cov)
		                                        if nameInd[i] in dicoCOV:
		                                            #print(nameInd[i])
		                                            #print(dicoCOV[nameInd[i]])
		                                            cutoffmini95CI = dicoCOV[nameInd[i]][2] # min cov quantile 5% 
		                                            cutoffmax95CI = dicoCOV[nameInd[i]][8] # max cov quantile 95%
		                                            #iprint(nameInd[i],cutoffmini95CI,cutoffmax95CI)
		                                            if int(cutoffmini95CI) > int(minCoverageThreshold):
		                                                minCoverageThreshold = int(cutoffmini95CI) # adjust the cutoff to min95CI
		                                            if int(cutoffmax95CI) < int(maxCoverageThreshold):
		                                                maxCoverageThreshold = int(cutoffmax95CI)
		                                            #print(nameInd[i],cutoffmini95CI,cutoffmax95CI,minCoverageThreshold,maxCoverageThreshold,minCoverageThresholdIni,maxCoverageThresholdIni)
							if cov >= int(minCoverageThreshold) and cov <= int(maxCoverageThreshold): # [;] rather than ];[
		                                                #myline="TRUE add=",str(REF),"cov=" + str(cov) + " within [" + str(minCoverageThreshold) + "-" + str(maxCoverageThreshold) + "]"
		                                                #print(myline)
								if ploidy[nameInd[numberindividual]]==2:
									#print("here")
									site.append(REF)
									site.append(REF)
								elif ploidy[nameInd[numberindividual]]==3:
									site.append(REF)
									site.append(REF)
									site.append(REF)
								elif ploidy[nameInd[numberindividual]]==4:
									site.append(REF)
									site.append(REF)
									site.append(REF)
									site.append(REF)
								else:
									break
							else:
		                                                #myline="FALSE add=N cov=" + str(cov) + " not within [" + str(minCoverageThreshold) + "-" + str(maxCoverageThreshold) + "]"
		                                                #print(myline)
								if ploidy[nameInd[numberindividual]]==2:
									#print("here")
									site.append("N")
									site.append("N")
								elif ploidy[nameInd[numberindividual]]==3:
									site.append("N")
									site.append("N")
									site.append("N")
								elif ploidy[nameInd[numberindividual]]==4:
									site.append("N")
									site.append("N")
									site.append("N")
									site.append("N")
								else:
									break
						numberindividual+=1
			else:   ## SNP SITE
				#print(len(sites)/NumberOfInd," ALT=",ALT)
				if arrline[6] != "PASS" or ALT == ".": # exclude lowQ SNP
					#print(len(sites)/NumberOfInd," ALT=",ALT)
					for i in range(0, len(nameSeq)): # already take care of the ploidy since it is for all nameSeq
						site.append("N")
						
					#print("lowQ  " + line)
				else:
					#print("PASS  " + line)
					gtPos = findPosInFORMATStr(arrline[8], "GT") #FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
					#gqPos = findPosInFORMATStr(arrline[8], "GQ") #FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">
					covPos = findPosInFORMATStr(arrline[8], "DP")
					numberindividual=0
					for i in range(0, len(nameInd)):
						ind = arrline[i+9]
						arrInd = ind.split(":")
						minCoverageThreshold = minCoverageThresholdIni
						maxCoverageThreshold = maxCoverageThresholdIni
						#print(nameInd[i])
						#print(dicoCOV)
						#print(numberindividual)
						if nameInd[i] in dicoCOV:
							#print(numberindividual,"=",nameInd[i]," found in dicoCOV")
                                               		cutoffmini95CI = dicoCOV[nameInd[i]][2]
                                                	cutoffmax95CI = dicoCOV[nameInd[i]][8]
							#print(minCoverageThreshold,maxCoverageThreshold)
                                                	if int(cutoffmini95CI) > int(minCoverageThreshold):
                                                    		minCoverageThreshold = int(cutoffmini95CI) # adjust the cutoff to min95CI
								#print("minCoverageThreshold changed min to",minCoverageThreshold)
                                                	if int(cutoffmax95CI) < int(maxCoverageThreshold):
                                                   		maxCoverageThreshold = int(cutoffmax95CI)
								#print("maxCoverageThreshold changed maximum coverage to ",maxCoverageThreshold)
                                                #print(dicoCOV[nameInd[i]])
                                                #print(minCoverageThreshold,maxCoverageThreshold,cutoffmini95CI,cutoffmax95CI)                                                                                                                                                               
						if len(arrInd) < int(covPos+1): # if individual doesn't have GQ field
							if ploidy[nameInd[numberindividual]]==2:
								#print("here <---------------------------- ploidy=",ploidy[nameInd[numberindividual]]),
								site.append("N")
								site.append("N")
							elif ploidy[nameInd[numberindividual]]==3:
								#print("here <----------------------------ploidy=",ploidy[nameInd[numberindividual]])
								site.append("N")
								site.append("N")
								site.append("N")
							elif ploidy[nameInd[numberindividual]]==4:
								#print("here <----------------------------ploidy=",ploidy[nameInd[numberindividual]])
								site.append("N")
								site.append("N")
								site.append("N")
								site.append("N")
							else:
								break
							#continue
							
						else:
							cov = arrInd[covPos]
							#print("coverage=",cov)
							#print("maxcoveragethreshold=",maxCoverageThreshold)
							if cov != "." and int(cov) <= int(maxCoverageThreshold) and int(cov) >= int(minCoverageThreshold):
								#print(nameInd[numberindividual]," coverage=",cov,"inf a",str(maxCoverageThreshold))
								arrInd = ind.split(":")
								#gt = arrInd[gtPos].split("/")
								#print(arrInd[gtPos])
								gt = re.split(r'/|\|', arrInd[gtPos])
								#print(gt)
								gt=shuffled(gt)
								#print(gt)
								if ploidy[nameInd[numberindividual]]==2:
									#print(str(nameInd[numberindividual]),"expected ploidy=",str(ploidy[nameInd[numberindividual]]),"=2?")
									if gt[0] == "0":
										site.append(REF)
									elif gt[0] == "1":
										site.append(ALT)
									else:
										site.append("N")
									if gt[1] == "0":
										site.append(REF)
									elif gt[1] == "1":
										site.append(ALT)
									else:
										site.append("N")
								elif ploidy[nameInd[numberindividual]]==3:
									#print(str(nameInd[numberindividual]),"expected ploidy=",str(ploidy[nameInd[numberindividual]]),"=3?")
									if gt[0] == "0":
										site.append(REF)
									elif gt[0] == "1":
										site.append(ALT)
									else:
										site.append("N")
									if gt[1] == "0":
										site.append(REF)
									elif gt[1] == "1":
										site.append(ALT)
									else:
										site.append("N")
									if gt[2] == "0":
										site.append(REF)
									elif gt[2] == "1":
										site.append(ALT)
									else:
										site.append("N")
								elif ploidy[nameInd[numberindividual]]==4:
									#print(str(nameInd[numberindividual]),"expected ploidy=",str(ploidy[nameInd[numberindividual]]),"=4?")
									if gt[0] == "0":
										site.append(REF)
									elif gt[0] == "1":
										site.append(ALT)
									else:
										site.append("N")
									if gt[1] == "0":
										site.append(REF)
									elif gt[1] == "1":
										site.append(ALT)
									else:
										site.append("N")
									if gt[2] == "0":
										site.append(REF)
									elif gt[2] == "1":
										site.append(ALT)
									else:
										site.append("N")
									if gt[3] == "0":
										site.append(REF)
									elif gt[3] == "1":
										site.append(ALT)
									else:
										site.append("N")
							else:
								#print(nameInd[numberindividual]," coverage=",cov,"sup a",str(maxCoverageThreshold))
								if ploidy[nameInd[numberindividual]]==2:
									#print("here")
									site.append("N")
									site.append("N")
								elif ploidy[nameInd[numberindividual]]==3:
									site.append("N")
									site.append("N")
									site.append("N")
								elif ploidy[nameInd[numberindividual]]==4:
									site.append("N")
									site.append("N")
									site.append("N")
									site.append("N")
								else:
									break
						numberindividual+=1
						#print("numberofindividual=",numberindividual)
		

		if len(sites) % 4.0 != 0.0 :
			#print("Site "+str(len(site)))
			#print("Sites "+str(len(sites)))
			print(str(pos))
			#sys.exit("NOOOO :-)")

		#sites.append(site)
		sites = sites + site
		oldpos = pos
		

print("Print  " + namefile)
fasta = open(namefile+".fst", "w")
sites = np.ravel(sites)
print(len(sites))
print(sites)
sequenceDNA = np.reshape(sites, newshape=(len(sites)/NumberOfInd, NumberOfInd))
for i, seqName in enumerate(nameSeq):
		print(">"+seqName, file=fasta)
		print(''.join(map(str, sequenceDNA[1:,i])), file=fasta)


#print "Number of SNP "+str(count)

fp.close()

