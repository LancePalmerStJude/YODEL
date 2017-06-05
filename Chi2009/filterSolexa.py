'''
Created on Jul 18, 2016

@author: lpalmer1
'''

'''

#band=110_20
band=130_50

for i in "A" "B" "C" "D" "E"
do

bsub -P peakcaller -q compbio_priority  -R "span[hosts=1] rusage[mem=8000]" \
python filterSolexa.py Brain${i}_${band}_fastq.txt  Brain${i}_${band}.filtered.fasta 
done



module load fastx_toolkit/0.0.13
module load cutadapt/1.2.1

for i in  "A" "B" "C" "D" "E"
do
bsub -P peakcaller -q compbio_priority  -R "span[hosts=1] rusage[mem=8000]" \
"cutadapt -a GTGTCAGTCACTTCCAGCGGTCGTATGCC -n 2 -m 16  -f fasta Brain${i}_${band}.filtered.fasta > Brain${i}_${band}.cutadapt.fasta " 
done



for i in  "A" "B" "C" "D" "E"
do
bsub -P peakcaller -q compbio_priority  -R "span[hosts=1] rusage[mem=8000]" \
fastx_collapser -i Brain${i}_${band}.cutadapt.fasta -o Brain${i}_${band}.collapsed.fasta
done



###for nonfiltered


####for not filtering
for i in "A" "B" "C" "D" "E"
do

bsub -P peakcaller -q compbio_priority  -R "span[hosts=1] rusage[mem=8000]" \
python filterSolexa.py Brain${i}_${band}_fastq.txt  Brain${i}_${band}.notfiltered.fasta 
done




for i in  "A" "B" "C" "D" "E"
do
bsub -P peakcaller -q compbio_priority  -R "span[hosts=1] rusage[mem=8000]" \
fastx_collapser -i Brain${i}_${band}.notfiltered.fasta -o Brain${i}_${band}.collapsednotfiltered.fasta
done

cat *collapsednotfiltered* > Brain_All_${band}.final.fasta
bsub -P HEMATO -q compbio_priority   -R "span[hosts=1] rusage[mem=32000]"  -oo logs/align.all.stdout -eo logs/align.all.stderr ~/programs/STAR/STAR-STAR_2.4.0k/bin/Linux_x86_64_static/STAR --genomeDir /nfs_exports/genomes/1/Mus_musculus/Mm9/STAR_small --genomeLoad NoSharedMemory  --readFilesIn  Brain_All_${band}.final.fasta  --outFileNamePrefix STARALL/Brain_All_${band}.star --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate  --quantMode TranscriptomeSAM  --chimSegmentMin 10 --chimJunctionOverhangMin 10 

mkdir COVERAGE
mkdir BED
mkdir CLUSTERS


samtools view -u  Brain_All_130_50.starAligned.sortedByCoord.out.bam | bamToBed -split -i stdin  > Brain_All_${band}.bed
LC_COLLATE=C sort  -k1,1 -k2,2n Brain_All_${band}.bed > Brain_All_${band}.sorted.bed
bedtools cluster -i Brain_All_${band}.sorted.bed -s  > Brain_All_${band}.clusters.bed
bedtools genomecov -split -bg  -i Brain_All_${band}.sorted.bed -g ~/data/sizes/mm9.sizes > Brain_All_${band}.sorted.cov
bedGraphToBigWig  Brain_All_${band}.sorted.cov  mm9.sizes  Brain_All_${band}.sorted.bw 


'''

###




import os 
import sys


inFile =sys.argv[1]
outFile=sys.argv[2]


IN=open(inFile)
OUT=open(outFile,'w')


while True:
    
    header=IN.readline()
    if not header:
        break
    seq=IN.readline()
    seq=seq.rstrip("\n")
    seq=seq.rstrip("\r")
    x=IN.readline()
    qual=IN.readline()
    qual=qual.rstrip("\n")
    qual=qual.rstrip("\r")
    qs=qual.split(" ")
    
    #seqs=seq.split()

    goodQ=1
    
    for i in range(len(qs)):
        q=qs[i]
        if int(q)<10:
            goodQ=0
            break
    if goodQ==0:
        continue
    OUT.write(">"+header[1:])
    OUT.write(seq+"\n")
    
    