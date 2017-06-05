'''
Created on Jul 22, 2016

@author: lpalmer1
'''
'''
Created on Jul 6, 2015

@author: lpalmer1
'''
#0 none
#5 single mismatch
#8 offset 
#11 6mer
#12 7mer-A1
#13 7mer-m8
#14 8mer

import sys
import os

import re
from Bio import SeqIO
from string import upper





seedFile=sys.argv[1]
clusterBase=sys.argv[2]
outBase=sys.argv[4]
outDir=sys.argv[3]



seedRegionsPosFile=os.path.join(clusterBase+".pos.fa")
seedRegionsNegFile=os.path.join(clusterBase+".neg.fa")

files=(['+',seedRegionsPosFile],['-',seedRegionsNegFile])


OUTBED=open(os.path.join(outDir,outBase+".seeds.bed"),'w')
OUTBEDPOS=open(os.path.join(outDir,outBase+".seeds.pos.bed"),'w')
OUTBEDNEG=open(os.path.join(outDir,outBase+".seeds.neg.bed"),'w')



SEEDS=open(seedFile)
header=SEEDS.readline()
seeds=[]
for line in SEEDS:
    line=line.rstrip("\n")
    line=line.rstrip("\r")
    if line[0]=="#":
        continue
    seedID,seedName,seedMatch,seed=line.split("\t")
    seeds.append([seedID,seedMatch])



for f in files:
    strand,cfile=f
    print(cfile)
    IN=open(cfile)
    for record in SeqIO.parse(IN, "fasta") :

        seq=record.seq
        sid=record.id
        chrm,se=sid.split(":")
        start,end=se.split("-")

        tmpSeq=upper(str(seq))
        #print(str(seq),strand)
        for sRef in seeds:
            mirName,seed=sRef
            mer6=seed[1:]
            #print(seed,tmpSeq)
             
            mer6PlusMatch=0
            for m in re.finditer(mer6, tmpSeq):
                matchClass=11
                mer6PlusMatch=1
                mStart=m.start()
                realStart=mStart
                length=6
                if mStart+7<=len(tmpSeq):
                    if tmpSeq[mStart+6]=="A":
                        matchClass=12
                        length=7
                if mStart>0:
                    if tmpSeq[mStart-1]==seed[0]:
                        realStart=mStart-1
                        if matchClass==12:
                            matchClass=14
                            length=8
                        else:
                            matchClass=13
                            length=7
                        
                
                
                startPos=int(start)+realStart
                
                
                if strand=='-':
                    startPos=int(start)+len(seq)-realStart-length
                endPos=startPos+length
                OUTBED.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
                if strand=="+":
                    OUTBEDPOS.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
                elif strand=="-":
                    OUTBEDNEG.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
            continue
            offSetMatch=0
            if mer6PlusMatch==0:
                mer6Offset=seed[0:6]
                for m in re.finditer(mer6Offset, tmpSeq):
                    offSetMatch=1
                    mStart=m.start()
                    realStart=mStart
                    startPos=int(start)+realStart
                    length=6
                    matchClass=8
                    if strand=='-':
                        startPos=int(start)+len(seq)-realStart-length
                    endPos=startPos+length
                    OUTBED.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
                    if strand=="+":
                        OUTBEDPOS.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
                    elif strand=="-":
                        OUTBEDNEG.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
            
            if offSetMatch==0 and mer6PlusMatch==0:
                if (mirName=="miR451a" or mirName == "miR1443p" or mirName=="miR1445p"):
                    for i in range(len(tmpSeq)-5):
                        sub=tmpSeq[i:i+6]
                        mismatches=0
                        mismatchPos=0
                        for j in range(6):
                            x=sub[j]
                            y=mer6[j]
                            if x!=y:
                                mismatches+=1
                                mismatchPos=6-j+1
                        if mismatches==1:
                            
                            realStart=i
                            startPos=int(start)+realStart
                            length=6
                            matchClass=mismatchPos
                            if strand=='-':
                                startPos=int(start)+len(seq)-realStart-length
                            endPos=startPos+length
                            OUTBED.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
                            if strand=="+":
                                OUTBEDPOS.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
                            elif strand=="-":
                                OUTBEDNEG.write("\t".join([chrm,str(startPos),str(endPos),mirName,str(matchClass),strand])+"\n")
                    
            
    IN.close()

