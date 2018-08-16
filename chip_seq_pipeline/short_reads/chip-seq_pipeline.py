#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob

bowtie_index = "/data/database/GRCh38/GENCODE/bowtie_index/GRCh38"
"""
# read length: 30~50bp
bowtie -q -m 1 -v 3 --sam --best --strata -p 10 /data/database/GRCh38/GENCODE/bowtie_index/GRCh38 input_mettl3.fq > input_mettl3.sam
bowtie -q -m 1 -v 3 --sam --best --strata -p 10 /data/database/GRCh38/GENCODE/bowtie_index/GRCh38 ip_mettl3.fq > ip_mettl3.sam

#
# ip
nohup samtools view -@ 10 -Sbh -o stat1_ip.bam stat1_ip.sam &
nohup samtools sort -@ 20 stat1_input.bam stat1_input_sorted &
rm unsorted_ip.bam
nohup samtools view -@ 20 -F4 -q 20 -b stat1_ip_sorted.bam -o stat1_ip_filtered.bam &
#
nohup java -jar /home/xiaoshan/bin/picard.jar MarkDuplicates I=stat1_ip_filtered.bam O=stat1_ip_unique.bam REMOVE_DUPLICATES=true CREATE_INDEX=true M=stat1_ip_rmdup.log &

#
# bedtools genomecov -ibam ip_mettl3_unique.bam -bga -g /data/database/hg38/hg38.chrom.sizes > ip_genomecov.bdg

macs2 callpeak -B -t ip_mettl3_unique.bam -c input_mettl3_unique.bam --outdir ../macs2_peak -q 0.05 --nomodel -n mettl3
cut -f 1-3 sample.narrowPeak > sample.peak.bed

#
for sample in input_mettl3_sorted ip_mettl3_sorted;do
    echo -en $sample"\t"
    raw = $(samtools view ${sample}.bam | wc -l)
    bamToBed -i ${sample}.bam | awk -vRAW = $raw ' {coordinates=$1":"$2"-"$3;
    total++;count[coordinates]++}END{for(coordinates in count){if(!max||count[coordinates]>max){max=count[coordinates];
    maxCoor=coordinates};if(count[coordinates]==1){unique++}};print
    RAW,total,total*100/RAW,unique,unique*100/total,maxCoor,count[maxCoor],count[maxCoor]*100/total}'
    samtools view -f 0x0004 ${sample}.bam | awk '{read=$10;total++;count[read]++}END{print "Total_non-mapped_reads",total;
    for(read in count){print read,count[read]+0}}' | sort -k2,2nr | head -11

done

"""


# short read < 50bp
def process_each_file(fq):
    result_sam = os.path.join(os.path.split(fq)[0], "result.sam")
    os.system("bowtie -q -m 1 -v 3 --sam --best --strata -p 10 %s %s > %s" % (bowtie_index, fq, result_sam))
    raw_bam = os.path.join(os.path.split(fq)[0], "result.sam")
    os.system("samtools view -@ 35 -Sbh -o unsorted_ip.bam ip_mettl3.sam")


def process_pipe_line():
    pass