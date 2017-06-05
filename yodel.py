#!/usr/bin/env python
description='''This program performs peak calling on clustered reads from clusterBed (bedTools).  
Peaks can be called from multiple sample at once, but samples will need to be
identified in readID from the input file.  The name should take the following
format sampleID:readID.  sampleID can contain any character except except white 
space and colon.  ReadID can be any character except white space.  A sample
file should also be given.  A header line is required and will ignored.  The
rest of the lines should contain sampleID a phenoType (which is ignored
currently) and finally 1 or 0 to indicate whether to use in peak calling (1) or
not (0).  These values should be tab delimited. Other parameters determine how
algorithm functions.
'''
__author__ = "Lance Palmer"
__copyright__ = "TBD"
__license__ = "TBD"
__version__ = "1.0"
__email__ = "Lance.Palmer@stjude.org"


import sys,getopt
from  os import path,makedirs
from operator import itemgetter
from decimal import Decimal


'''
--mph=5 --inFile=C:\projects\hitsClip\PeakCaller\allSamples.fullCollapsed.clusters.bed --outDir=C:\projects\hitsClip\PeakCaller --dt=1.5 --pdb=1 --prefix=MoreSensTest  --sampleList=C:\projects\hitsClip\PeakCaller\sampleList.txt
'''

####Default values####
#minmimum number of reads in a peak
minPeakHeight=3
#See formula for peak calling
peakDipBuffer=1
dipTolerance=1.5
#size
binSize=16
#input file from clusterBed.  Feature name in bed file should contain sample id first followed by ':' then by sequence id
clusterFile=''
#output directory
outDir=''
#list of samples.  After header each line should be tab delimited containing sampleName\tPhenotype\tUsage (1 to use for peak calling, 0 for not)
sampleListFile=''
#--sampleList=C:\projects\hitsClip\PeakCaller\sampleList.txt
#prefix for output files
outPrefix="peaks"
#wether sample list is present. if it is not, peaks made on all reads
gotSampleList=0
usage='''peakCaller -i <clusterFile.bed> -o <output directory> -s <sampleList> [--mph --pdb --dt --bs --prefix]
--mph        minPeakCount  Default=6
--pdb        peakDipBuffer  Default=1
--dt         dipTolerance   Default=1.5
--bs         binSize        Default=16
--prefix     output prefix  Default=peaks
''' +description
###############Parse parameters
try:
    opts, args = getopt.getopt(sys.argv[1:],'i:o:s:',["sampleList=","mph=","pdb=","inFile=","outDir=","dt=","bs=","prefix="])
except getopt.GetoptError as err:
    print err
    sys.exit()
try:
    for opt, arg in opts:  
        if opt=="--mph":
            minPeakHeight=int(arg)
        elif opt=="--pdb":
            peakDipBuffer=int(arg)
        elif opt=="--dt":
            dipTolerance=Decimal(arg)
        elif opt=="--bs":
            binSize=int(arg)
        elif opt in ("-i","--inFile"):
            clusterFile=arg
        elif opt in ("-o","--outDir"):
            outDir=arg     
        elif opt in ("-s","--sampleList"):
            sampleListFile=arg
        elif opt=="--prefix":
            outPrefix=arg
except StandardError as err:
    print(err)
    print(usage)
    sys.exit(2)
            
if len(clusterFile)==0:
    print("inFile not defined")
    print(usage)
    sys.exit(2)
if not path.isfile(clusterFile):
    print("inFile does not exist")
    print(usage)
    sys.exit(2)
if len(outDir)==0:
    print("outDir not defined")
    print(usage)
    sys.exit(2)
try:
    if not path.exists(outDir):
        makedirs(outDir)   
except:
    print("Could not make directory:"+outDir)
if len(sampleListFile)>0:
    gotSampleList=1
    if not path.isfile(sampleListFile):
        print("sampleListFile does not exist")
        print(usage)
        sys.exit(2)

    

def printParameters():
    print("Parameters:")
    print("inFile\t"+clusterFile)
    print("outDir\t"+outDir)
    print("sampleListFile\t"+sampleListFile)
    print("minPeakHeight\t"+str(minPeakHeight))
    print("peakDipBuffer\t"+str(peakDipBuffer))
    print("dipTolerance\t"+str(dipTolerance))
    print("binSize\t"+str(binSize))    
    print("prefix\t"+outPrefix) 
printParameters()
#End parsing parameters############################################
binCount=0
sampleList=[]
samples={}
phenos={}
#parse sample file

if gotSampleList==1:
    SAMPLE=open(sampleListFile)
    sampleHeader=SAMPLE.readline()
    #3 columns for sample file  SampleName,pheno,wether to use in peak calling (1 for yes , 0 for no)
    # example : 1WT    WT    1
    for line in SAMPLE:
        line=line.rstrip("\n")
        line=line.rstrip("\r")
        sid,pheno,use=line.split("\t")
        samples[sid]=[pheno,use,0] # index 2 is a count
        sampleList.append(sid)
        if pheno not in phenos:
            phenos[pheno]=[]
        phenos[pheno].append(sid)
else:
    #use all samples 
    sampleList.append("All")
    samples["All"]=["All","1",0]

#######Output handles
outFile = path.join(outDir, outPrefix + ".bins.cov.bed")
OUT = open(outFile, 'w')
outFileFull = path.join(outDir, outPrefix + ".bins.cov.full.bed")
OUTFULL = open(outFileFull, 'w')
outBinBounds = path.join(outDir, outPrefix + ".binBounds.txt")
BINBOUNDS = open(outBinBounds, 'w')
binBoundsData = ("binName", "chr", "binStart", "binEnd", "topHeight", "topPos", "start25", "end25", "start50", "end50", "start75", "end75", "peakHeight", "readCounts25","readCounts50")
BINBOUNDS.write("\t".join(binBoundsData) + "\n")
if gotSampleList==1:
    outTableRC25 = path.join(outDir, outPrefix + ".binCountsRC25.txt")
    OUTTABLERC25 = open(outTableRC25, 'w')
    outTableRC50 = path.join(outDir, outPrefix + ".binCountsRC50.txt")
    OUTTABLERC50 = open(outTableRC50, 'w')
    outTablePH = path.join(outDir, outPrefix + ".binCountsPH.txt")
    OUTTABLEPH = open(outTablePH, 'w')
    sHeader="Bin"
    for sid in sampleList:
        sHeader+="\t"+sid
    OUTTABLERC25.write(sHeader+"\n")   
    OUTTABLERC50.write(sHeader+"\n")   
    OUTTABLEPH.write(sHeader+"\n") 


#input is start and end of interval 1 and start and end of interval 2
#output=positive number is amount of overlap, negative number is distance
def overlap(s1, e1, s2, e2):
    if s1 > e1:
        tmp = s1
        s1 = e1
        e1 = tmp
    if s2 > e2:
        tmp = s2
        s2 = e2
        e2 = tmp   
    dist = 0
    if s2 < s1:
        tmps = s2
        s2 = s1
        s1 = tmps
        tmpe = e2
        e2 = e1
        e1 = tmpe
    if e2 >= e1:
        dist = e1 - s2
        if dist > 0:
            dist = dist + 1
    else:
        dist = e2 - s2 + 1
    return dist

def processPeaks(tmpReads):
    global binCount
    global sampleCounts
    #baseCoverage=at each base number of reads for samples being used fore peak calling
    #baseCoverageBySample coverage at each base for all samples
    #goodReads number of usable reads in regions, only used to compare to minpeakheight
    #readInfo contains start and end of each read (1 based)
    #readID2Sid links read id to sample
    #cluster varoane not used... its the cluster number form input bed file
    #chrm and strand are the chromosome and strand of the cluster of reads
    baseCoverage, baseCoverageBySample, goodReads, readInfo, readID2Sid, cluster, chrm, strand = getCoverage(tmpReads)    
    if goodReads < minPeakHeight:
        return
    #get bounds of cluster
    sortByBase = sorted(baseCoverage.items(), key=itemgetter(0))
    lowPos = sortByBase[0][0]
    hiPos = sortByBase[-1][0]
    #as reads removed when peaks are built, tmpCoverage will change
    tmpCoverage = baseCoverage
    #hold final peak info
    peaks = []
    #print(tmpCoverage)
    while(len(tmpCoverage) > 0):
        #to get position with highest coverage
        sortedCov = sorted(tmpCoverage.items(), key=itemgetter(1), reverse=True)
        #position of highest coverage for peak
        topPos = sortedCov[0][0]
        #actual number of reads at that position
        topHeight = sortedCov[0][1]
        if topHeight < 1:
            break
        if topHeight<minPeakHeight:
            break
        #if there are multiple maximum for  apeak record all those positions, to deal with later
        maxHeightMulti=[]
        maxHeightMulti.append(topPos)
        #will recored lowest position going downstream of highest point
        minDown = topHeight
        #record the position of lowest point
        minDownPos = topPos
        #denote peak boundaries
        peakBegin = lowPos
        #go from heighest peak position upstream tell end of peak found.  lowPos is beginning of cluster
        for i in range(topPos - 1, lowPos - 2, -1):
            height = 0
            if i in tmpCoverage:
                height = tmpCoverage[i]
            #calculate if peak should end based on height of peak at that position
            if height  >= (minDown + peakDipBuffer) * dipTolerance or height == 0:
                peakBegin = minDownPos
                break
            #keep track of lowest position
            if height <= minDown:
                minDown = height
                minDownPos = i
            #if multiple top positions
            if height==topHeight:
                maxHeightMulti.append(i) 
        #now for downstream of peak (up in position)  
        minUp = topHeight
        minUpPos = topPos
        peakEnd = hiPos
        # see for loop above for comments
        for i in range(topPos + 1, hiPos + 2):
            height = 0
            if i in tmpCoverage:
                height = tmpCoverage[i]
            if height  >= (minUp + peakDipBuffer) * dipTolerance or height == 0:
                peakEnd = minUpPos
                break
            if height <= minUp:
                minUp = height
                minUpPos = i  
            if height==topHeight:
                maxHeightMulti.append(i)    
        #these will represent coordinates at 25%,50% and 75% of peak height
        up75 = -1
        up50 = -1
        up25 = -1 
        down75 = -1
        down50 = -1
        down25 = -1              
        #getting coordinates of 25%,50% and 75% of peak height
        for i in range(topPos, peakEnd + 1):
            height = tmpCoverage[i]
            if height >= topHeight * .75:
                up75 = i
            if height >= topHeight * .50:
                up50 = i
            if height >= topHeight * .25:
                up25 = i                
        for i in range(topPos, peakBegin - 1, -1):
            height = tmpCoverage[i]
            if height >= topHeight * .75:
                down75 = i
            if height >= topHeight * .50:
                down50 = i
            if height >= topHeight * .25:
                down25 = i         
        #if multiple high points, take one in middle 
        sortTopHeights=sorted(maxHeightMulti)
        midTop=sortTopHeights[int(len(sortTopHeights)/2)]    
        #storing peak info  
        peaks.append([peakBegin, peakEnd, topHeight,midTop, down25, up25, down50, up50, down75, up75])
        #set all coverag in peak to 0
        for i in range(peakBegin, peakEnd + 1):
            tmpCoverage[i] = 0
    #loop through all peaks found
    for p in peaks:
        peakBegin, peakEnd, topHeight, topPos,down25, up25, down50, up50, down75, up75 = p
        #if peaksize is less then binSize, do not keep
        if up25 - down25 + 1 < binSize:
            continue
        #bin will be named based on 25% of peak height marks
        binName = chrm + ":" + str(down25) + ":" + str(up25) + ":" + strand
        #stores peak height
        topCountsPH = {} 
        #stores read count for peak 
        topCountsRC25 = {}
        topCountsRC50 = {}
        #go samplecoverage,  Calculate coverage per sample.  Use all samples, not just for calling peaks
        for sid in sampleList:
            topHeightSid = 0
            topCountsRC25[sid] = 0
            topCountsRC50[sid] = 0
            for i in range(down25, up25 + 1):
                #count number of positions covered by read that overlap peak
                if i in baseCoverageBySample[sid]:              
                    count = baseCoverageBySample[sid][i]
                    #record highest count
                    if count > topHeightSid:
                        topHeightSid = count
            #highest peak for sample within peak
            topCountsPH[sid] = topHeightSid        
        #go through reads, calculate number of reads that overlap 25 and 50% peaks by at least bin size
        for readID in readInfo:
            rstart, rend = readInfo[readID]            
            ov25 = overlap(rstart, rend, down25, up25)
            ov50 = overlap(rstart, rend, down50, up50)
            if  ov25 >= binSize:
                sid = readID2Sid[readID]
                #number of reads overlapping
                topCountsRC25[sid] += 1
            if  ov50 >= binSize:
                sid = readID2Sid[readID]
                #number of reads overlapping
                topCountsRC50[sid] += 1 
        ###preparing output
        outLinePH = binName
        outLineRC25 = binName
        outLineRC50 = binName
        totalCountsPH = 0
        totalCountsRC25 = 0
        totalCountsRC50 = 0
        #per sample outputs
        for sid in sampleList:
            outLinePH = outLinePH + "\t" + str(topCountsPH[sid])
            outLineRC25 = outLineRC25 + "\t" + str(topCountsRC25[sid])
            outLineRC50 = outLineRC50 + "\t" + str(topCountsRC50[sid])
            if samples[sid][1]!="1":
                continue
            totalCountsRC25 += topCountsRC25[sid]
            totalCountsRC50 += topCountsRC25[sid]            
            totalCountsPH += topCountsPH[sid] 
        if topHeight < minPeakHeight:
            return
        ##next 3 files, stats for differential peaks
        if gotSampleList==1:
            OUTTABLEPH.write(outLinePH + "\n")
            OUTTABLERC25.write(outLineRC25 + "\n")
            OUTTABLERC50.write(outLineRC50 + "\n")
        binCount += 1     
        # basic bin stats
        OUT.write("\t".join([chrm, str(down25 - 1), str(up25), binName, str(totalCountsRC25), strand, str(down50 - 1), str(up50)]) + "\n")
        OUTFULL.write("\t".join([chrm, str(peakBegin), str(peakEnd), strand, str(totalCountsRC25), strand,str(down25 - 1), str(up25)]) + "\n")
        #bin info
        binBoundsData = (binName, chrm, peakBegin, peakEnd, topHeight, topPos, down25, up25, down50, up50, down75, up75, str(totalCountsPH), str(totalCountsRC25),str(totalCountsRC25))
        BINBOUNDS.write("\t".join(map(str, binBoundsData)) + "\n") 

#for each cluster, get coverage at each position for all samples and samples for peak calling
def getCoverage(tmpReads):
    #start end for each read
    readInfo = {}
    #link read to sample id
    readID2Sid = {}
    #base coverage at each position for peak calling samples
    baseCoverage = {}
    #count number of reads used for peaks (if reads only from controls, no need to do further peak calling
    goodReads = 0
    chrm = ''
    strand = ''
    cluster = -1
    baseCoverageBySample = {}
    for sid in sampleList:
        baseCoverageBySample[sid] = {}
    for r in tmpReads:
        sid='All'
        chrm, rstart, rend, readid, x, strand, cluster = r
        readInfo[readid] = [int(rstart), int(rend)]  # use true start
        if gotSampleList==1:
            rids = readid.split(":")  #sample id must be first
            sid = rids[0]
        readID2Sid[readid] = sid

        #go through coordinates of all reads and populate baseCoverageBySample
        for i in range(rstart, rend + 1):
            if i not in baseCoverageBySample[sid]:
                baseCoverageBySample[sid][i] = 0
            baseCoverageBySample[sid][i] += 1
        #only continue with samples to be used for base calling
        if samples[sid][1]!="1":
            continue
        #redo counting, all samples at once
        for i in range(rstart, rend + 1):
            if i not in baseCoverage:
                baseCoverage[i] = 0
            baseCoverage[i] += 1
            
        goodReads += 1
    return(baseCoverage, baseCoverageBySample, goodReads, readInfo, readID2Sid, cluster, chrm, strand)
    
CF = open(clusterFile)
tmpReads = []
lastCluster = None
#parse through output of 
for line in CF:
    line = line.rstrip("\n")
    line = line.rstrip("\r")
    chrm, start, end, readid, x, strand, cluster = line.split("\t")
    

    #should be form of sampleID:readID (read id can have : in it)
    rids = readid.split(":")
    sid = rids[0]
    #if read size is less than binsize, next
    if ((int(end) - int(start)) < int(binSize)):
        continue
    #first read will need to set cluster
    if lastCluster is None:
        lastCluster = cluster
    #if new cluster found
    if cluster != lastCluster:
        processPeaks(tmpReads)
        tmpReads = []
        lastCluster = cluster
    #store read info until new cluster found
    tmpReads.append([chrm, int(start) + 1, int(end), readid, x, strand, cluster])

processPeaks(tmpReads) #for that last cluster
print("Bins found",binCount)