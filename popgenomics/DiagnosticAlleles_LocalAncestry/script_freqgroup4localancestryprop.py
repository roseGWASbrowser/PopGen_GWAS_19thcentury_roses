#! /usr/bin/env python
from __future__ import print_function
import sys
import re

namevcf=sys.argv[1]
vcf = open(namevcf)
nameoutput=sys.argv[2]
outfreq = open(nameoutput,"w")
nameInd = []
PASSsite=0

for line in vcf:
	line = line.rstrip()
	if line[0:6] == "#CHROM" and len(nameInd) == 0:
		arrline = line.split()
		for i in arrline[9:]:
			nameInd.append(i)
		print(nameInd)
		myline=str("chromosome")+"\t"+str("position")+"\t"+str("REF")+"\t"+str("ALT")+"\t"+str("freqASIAcountref")+"\t"+str("freqEUROcountref")+"\t"+str("freqPHYBcountref")+"\t"+str("freqHYBTcountref")+"\t"+str("ASIAcountref")+"\t"+str("ASIAcountalt")+"\t"+str("EUROcountref")+"\t"+str("EUROcountalt") + "\t" + str("PHYBcountref") +"\t"+ str("PHYBcountalt") + "\t" + str("HYBTcountref")+"\t"+ str("HYBTcountalt") +"\n"
	elif line[0] != "#":
		arrline = line.split()
		#print(arrline)
		chromosome = arrline[0]
		pos = int(arrline[1])
		REF = str(arrline[3])
		ALT = str(arrline[4])
		if arrline[6] == "PASS":
			#print(line)
			if len(REF) > 1 or len(ALT) > 1 or ALT == "*" : 
				continue
			else:
				PASSsite+=1
				ASIAcountref=0
				ASIAcountalt=0
				ASIAcounttotal=0
				ASIAcountmissing=0
				EUROcountref=0
				EUROcountalt=0
				EUROcounttotal=0
				EUROcountmissing=0
				PHYBcountref=0
				PHYBcountalt=0
				PHYBcounttotal=0
				PHYBcountmissing=0
				HYBTcountref=0
				HYBTcountalt=0
				HYBTcounttotal=0
				HYBTcountmissing=0
				for i in range(0, len(nameInd)):
					#print(nameInd[i])
					ind = arrline[i+9]
					arrInd = ind.split(":")
					gt = re.split(r'/|\|', arrInd[0])
					#print(gt)
					#print(len(gt))
					if nameInd[i]== "175584" or nameInd[i]=="175585" or nameInd[i]=="OB2n" or nameInd[i]=="Rchimut" or nameInd[i]=="Rchisan" or nameInd[i]=="Rchispo" or nameInd[i]=="Rogig":
						if len(gt)==2:
							if gt[0] == "0":
								ASIAcountref+=1
								ASIAcounttotal+=1
							elif gt[0]=="1":
								ASIAcountalt+=1
								ASIAcounttotal+=1
							elif gt[0]==".":
								ASIAcountmissing+=1
								ASIAcounttotal+=1
							if gt[1] == "0":
								ASIAcountref+=1
								ASIAcounttotal+=1
							elif gt[1]=="1":
								ASIAcountalt+=1
								ASIAcounttotal+=1
							elif gt[1]==".":
								ASIAcountmissing+=1
								ASIAcounttotal+=1
						if len(gt)==4:
							if gt[0] == "0":
								ASIAcountref+=1
								ASIAcounttotal+=1
							elif gt[0]=="1":
								ASIAcountalt+=1
								ASIAcounttotal+=1
							elif gt[0]==".":
								ASIAcountmissing+=1
								ASIAcounttotal+=1
							if gt[1] == "0":
								ASIAcountref+=1
								ASIAcounttotal+=1
							elif gt[1]=="1":
								ASIAcountalt+=1
								ASIAcounttotal+=1
							elif gt[1]==".":
								ASIAcountmissing+=1
								ASIAcounttotal+=1
							if gt[2] == "0":
								ASIAcountref+=1
								ASIAcounttotal+=1
							elif gt[2]=="1":
								ASIAcountalt+=1
								ASIAcounttotal+=1
							elif gt[2]==".":
								ASIAcountmissing+=1
								ASIAcounttotal+=1
							if gt[3] == "0":
								ASIAcountref+=1
								ASIAcounttotal+=1
							elif gt[3]=="1":
								ASIAcountalt+=1
								ASIAcounttotal+=1
							elif gt[3]==".":
								ASIAcountmissing+=1
								ASIAcounttotal+=1
					if nameInd[i]== "175581" or nameInd[i]=="175583" or nameInd[i]=="177142" or nameInd[i]=="Rgaloff":
						if len(gt)==4:
							if gt[0] == "0":
								EUROcountref+=1
								EUROcounttotal+=1
							elif gt[0]=="1":
								EUROcountalt+=1
								EUROcounttotal+=1
							elif gt[0]==".":
								EUROcountmissing+=1
								EUROcounttotal+=1
							if gt[1] == "0":
								EUROcountref+=1
								EUROcounttotal+=1
							elif gt[1]=="1":
								EUROcountalt+=1
								EUROcounttotal+=1
							elif gt[1]==".":
								EUROcountmissing+=1
								EUROcounttotal+=1
							if gt[2] == "0":
								EUROcountref+=1
								EUROcounttotal+=1
							elif gt[2]=="1":
								EUROcountalt+=1
								EUROcounttotal+=1
							elif gt[2]==".":
								EUROcountmissing+=1
								EUROcounttotal+=1
							if gt[3] == "0":
								EUROcountref+=1
								EUROcounttotal+=1
							elif gt[3]=="1":
								EUROcountalt+=1
								EUROcounttotal+=1
							elif gt[3]==".":
								EUROcountmissing+=1
								EUROcounttotal+=1
					if nameInd[i]== "175586" or nameInd[i]=="175587" or nameInd[i]=="175623" or nameInd[i]=="175624":
						if len(gt)==4:
							if gt[0] == "0":
								PHYBcountref+=1
								PHYBcounttotal+=1
							elif gt[0]=="1":
								PHYBcountalt+=1
								PHYBcounttotal+=1
							elif gt[0]==".":
								PHYBcountmissing+=1
								PHYBcounttotal+=1
							if gt[1] == "0":
								PHYBcountref+=1
								PHYBcounttotal+=1
							elif gt[1]=="1":
								PHYBcountalt+=1
								PHYBcounttotal+=1
							elif gt[1]==".":
								PHYBcountmissing+=1
								PHYBcounttotal+=1
							if gt[2] == "0":
								PHYBcountref+=1
								PHYBcounttotal+=1
							elif gt[2]=="1":
								PHYBcountalt+=1
								PHYBcounttotal+=1
							elif gt[2]==".":
								PHYBcountmissing+=1
								PHYBcounttotal+=1
							if gt[3] == "0":
								PHYBcountref+=1
								PHYBcounttotal+=1
							elif gt[3]=="1":
								PHYBcountalt+=1
								PHYBcounttotal+=1
							elif gt[3]==".":
								PHYBcountmissing+=1
								PHYBcounttotal+=1
					if nameInd[i]== "175624" or nameInd[i]=="175627" or nameInd[i]=="175628" or nameInd[i]=="177143":
						if len(gt)==4:
							if gt[0] == "0":
								HYBTcountref+=1
								HYBTcounttotal+=1
							elif gt[0]=="1":
								HYBTcountalt+=1
								HYBTcounttotal+=1
							elif gt[0]==".":
								HYBTcountmissing+=1
								HYBTcounttotal+=1
							if gt[1] == "0":
								HYBTcountref+=1
								HYBTcounttotal+=1
							elif gt[1]=="1":
								HYBTcountalt+=1
								HYBTcounttotal+=1
							elif gt[1]==".":
								HYBTcountmissing+=1
								HYBTcounttotal+=1
							if gt[2] == "0":
								HYBTcountref+=1
								HYBTcounttotal+=1
							elif gt[2]=="1":
								HYBTcountalt+=1
								HYBTcounttotal+=1
							elif gt[2]==".":
								HYBTcountmissing+=1
								HYBTcounttotal+=1
							if gt[3] == "0":
								HYBTcountref+=1
								HYBTcounttotal+=1
							elif gt[3]=="1":
								HYBTcountalt+=1
								HYBTcounttotal+=1
							elif gt[3]==".":
								HYBTcountmissing+=1
								HYBTcounttotal+=1
				if ASIAcountmissing <= 4:
					freqASIAcountref=float(ASIAcountref)/(ASIAcountref+ASIAcountalt)
				else:
					freqASIAcountref="NA"
				if EUROcountmissing <= 4:
					freqEUROcountref=float(EUROcountref)/(EUROcountref+EUROcountalt)
				else:
					freqEUROcountref="NA"
				if PHYBcountmissing <= 4:
					freqPHYBcountref=float(PHYBcountref)/(PHYBcountref+PHYBcountalt)
				else:
					freqPHYBcountref="NA"
				if HYBTcountmissing <= 4:
					freqHYBTcountref=float(HYBTcountref)/(HYBTcountref+HYBTcountalt)
				else:
					freqHYBTcountref="NA"
				myline=str(chromosome)+"\t"+str(pos)+"\t"+str(REF)+"\t"+str(ALT)+"\t"+str(freqASIAcountref)+"\t"+str(freqEUROcountref)+"\t"+str(freqPHYBcountref)+"\t"+str(freqHYBTcountref)+"\t"+str(ASIAcountref)+"\t"+str(ASIAcountalt)+"\t"+str(EUROcountref)+"\t"+str(EUROcountalt) + "\t" + str(PHYBcountref) +"\t"+ str(PHYBcountalt) + "\t" + str(HYBTcountref)+"\t"+ str(HYBTcountalt) +"\n"
				#print(myline)
	else:
		print("unexpected situation in line:",line)
		break
	outfreq.write(myline)

outfreq.close()
vcf.close()				
				#print(PASSsite)
				#print(ASIAcountmissing," ", EUROcountmissing)
				#print(ASIAcountref)
				#print(ASIAcountalt)
				#print(ASIAcounttotal)
