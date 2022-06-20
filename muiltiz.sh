#!/bin/bash
set -e
set -u
set -o pipefail

# set the work path
work_path=$1

# make the tree file and test
mkdir -p ${work_path}/result
mkdir -p ${work_path}/result/00_tmp
mkdir -p ${work_path}/result/01_masked_genome
mkdir -p ${work_path}/result/02_2bit_chromesize
mkdir -p ${work_path}/result/03_pair_lastz
mkdir -p ${work_path}/result/03_filter_alignment


ls ${work_path}/genome/*.fa | while read id;
do
        readlink -f ${id}  >> ${work_path}/result/00_tmp/genomelist.txt
done

masked_genome=${work_path}/result/01_masked_genome

#step1:reapetmask,trf mask
echo trf begins
cat ${work_path}/result/00_tmp/genomelist.txt | while read id;
do
    base=$(basename ${id} .fa)
    echo Sample name is ${base}
    trf ${id} 2 7 7 80 10 50 500 -f -d -m -h
    mv ${base}.*.mask ${masked_genome}/${base}.mask.fa
    echo ${base} is ok
done

rm *.dat

ls ${work_path}/result/01_masked_genome/*.mask.fa | while read id;
do
        readlink -f ${id} >> ${work_path}/result/00_tmp/masked_genomelist.txt
done

twobit_chrsz=${work_path}/result/02_2bit_chromesize
#step2:transform to 2bit format and get chromesize
echo twobit_chrsz begins
cat ${work_path}/result/00_tmp/masked_genomelist.txt | while read id;
do
    base=$(basename ${id} .mask.fa)
    echo Sample name is ${base}
    faToTwoBit ${id} ${twobit_chrsz}/${base}.2bit
    twoBitInfo ${twobit_chrsz}/${base}.2bit stdout | sort -k2,2nr > ${twobit_chrsz}/${base}.chromsizes
    echo ${base} is ok
done


pair_align=${work_path}/result/03_pair_lastz
#step3:pairwise alignment,arabidopsis_thaliana as reference
echo pairwise alignment begins
cat ${work_path}/result/00_tmp/masked_genomelist.txt | while read id;
do
    base=$(basename ${id} .mask.fa)
    echo Sample name is ${base}
    lastz-1.04.00 ${masked_genome}/arabidopsis_thaliana.mask.fa[multiple] ${id}[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=axt > at_vs_${base}.axt
    echo ${base} is ok
done;



#step4:format transform axt to psl
nohup axtToPsl at_vs_at.chain.mask.axt arabidopsis_thaliana.chrsize arabidopsis_thaliana.chrsize at_vs_at.psl &
nohup axtToPsl at_vs_eug.chain.mask.axt arabidopsis_thaliana.chrsize eucalyptus_grandis.chrsize at_vs_eug.psl &
nohup axtToPsl at_vs_ors.chain.mask.axt arabidopsis_thaliana.chrsize oryza_sativa.chrsize at_vs_ors.psl &
nohup axtToPsl at_vs_pot.chain.mask.axt arabidopsis_thaliana.chrsize populus_trichocarpa.chrsize at_vs_pot.psl &
nohup axtToPsl at_vs_ra.chain.mask.axt arabidopsis_thaliana.chrsize rhizophora_apiculata.chrsize at_vs_ra.psl &
nohup axtToPsl at_vs_sa.chain.mask.axt arabidopsis_thaliana.chrsize sonneratia_alba.chrsize at_vs_sa.psl &
nohup axtToPsl at_vs_sei.chain.mask.axt arabidopsis_thaliana.chrsize se.indicum.chrsize at_vs_sei.psl &：
nohup axtToPsl at_vs_am.chain.mask.axt arabidopsis_thaliana.chrsize avicennia_marina.chrsize at_vs_am.psl &


#step5:axtChain
nohup axtChain -linearGap=loose -psl at_vs_ra.psl -faT arabidopsis_thaliana.fa -faQ rhizophora_apiculata.fa at_vs_ra.chain &
nohup axtChain -linearGap=loose -psl at_vs_sa.psl -faT arabidopsis_thaliana.fa -faQ sonneratia_alba.fa at_vs_sa.chain &
nohup axtChain -linearGap=loose -psl at_vs_pot.psl -faT arabidopsis_thaliana.fa -faQ populus_trichocarpa.fa at_vs_pot.chain &
nohup axtChain -linearGap=loose -psl at_vs_eug.psl -faT arabidopsis_thaliana.fa -faQ eucalyptus_grandis.fa at_vs_eug.chain &
nohup axtChain -linearGap=loose -psl at_vs_sei.psl -faT arabidopsis_thaliana.fa -faQ se.indicum.fa at_vs_sei.chain &
nohup axtChain -linearGap=loose -psl at_vs_ors.psl -faT arabidopsis_thaliana.fa -faQ oryza_sativa.fa at_vs_ors.chain &
nohup axtChain -linearGap=loose -psl at_vs_am.psl -faT arabidopsis_thaliana.fa -faQ avicennia_marina.fa at_vs_am.chain &
nohup axtChain -linearGap=loose -psl at_vs_at.psl -faT arabidopsis_thaliana.fa -faQ arabidopsis_thaliana.fa at_vs_at.chain &


#step6:Netting,chainNet
#a:chainMergeSort
for i in *.chain;
do
        chainMergeSort $i > $i.chain.sort
done;

#b:chainPreNet
nohup chainPreNet at_vs_at.chain.chain.sort arabidopsis_thaliana.chrsize arabidopsis_thaliana.chrsize at_vs_at.prenet &
nohup chainPreNet at_vs_eug.chain.chain.sort arabidopsis_thaliana.chrsize eucalyptus_grandis.chrsize at_vs_eug.prenet &
nohup chainPreNet at_vs_ors.chain.chain.sort arabidopsis_thaliana.chrsize oryza_sativa.chrsize at_vs_ors.prenet &
nohup chainPreNet at_vs_pot.chain.chain.sort arabidopsis_thaliana.chrsize populus_trichocarpa.chrsize at_vs_pot.prenet &
nohup chainPreNet at_vs_ra.chain.chain.sort arabidopsis_thaliana.chrsize rhizophora_apiculata.chrsize at_vs_ra.prenet &
nohup chainPreNet at_vs_sa.chain.chain.sort arabidopsis_thaliana.chrsize sonneratia_alba.chrsize at_vs_sa.prenet &
nohup chainPreNet at_vs_sei.chain.chain.sort arabidopsis_thaliana.chrsize se.indicum.chrsize at_vs_sei.prenet &
nohup chainPreNet at_vs_am.chain.chain.sort arabidopsis_thaliana.chrsize avicennia_marina.chrsize at_vs_am.prenet &

#c:chainNet
nohup chainNet at_vs_am.prenet arabidopsis_thaliana.chrsize avicennia_marina.chrsize atam.net am.net &
nohup chainNet at_vs_sa.prenet arabidopsis_thaliana.chrsize sonneratia_alba.chrsize atsa.net sa.net &
nohup chainNet at_vs_ra.prenet arabidopsis_thaliana.chrsize rhizophora_apiculata.chrsize atra.net ra.net &
nohup chainNet at_vs_ors.prenet arabidopsis_thaliana.chrsize oryza_sativa.chrsize ators.net ors.net &
nohup chainNet at_vs_pot.prenet arabidopsis_thaliana.chrsize populus_trichocarpa.chrsize atpot.net pot.net &
nohup chainNet at_vs_sei.prenet arabidopsis_thaliana.chrsize se.indicum.chrsize atsei.net sei.net &
nohup chainNet at_vs_eug.prenet arabidopsis_thaliana.chrsize eucalyptus_grandis.chrsize ateug.net eug.net &
nohup chainNet at_vs_at.prenet arabidopsis_thaliana.chrsize arabidopsis_thaliana.chrsize atat.net at.net &

#d:netSyntenic
for i in *.net; do netSyntenic $i $i.netsyn;done;

#step7:maffting
#a:netToAxt
nohup netToAxt atam.net.netsyn at_vs_am.prenet arabidopsis_thaliana.mask.2bit avicennia_marina.mask.2bit atam.net2axt &
nohup netToAxt atat.net.netsyn at_vs_at.prenet arabidopsis_thaliana.mask.2bit arabidopsis_thaliana.mask.2bit atat.net2axt &
nohup netToAxt ateug.net.netsyn at_vs_eug.prenet arabidopsis_thaliana.mask.2bit eucalyptus_grandis.mask.2bit ateug.net2axt &
nohup netToAxt ators.net.netsyn at_vs_ors.prenet arabidopsis_thaliana.mask.2bit oryza_sativa.mask.2bit ators.net2axt &
nohup netToAxt atpot.net.netsyn at_vs_pot.prenet arabidopsis_thaliana.mask.2bit populus_trichocarpa.mask.2bit atpot.net2axt &
nohup netToAxt atsei.net.netsyn at_vs_sei.prenet arabidopsis_thaliana.mask.2bit se.indicum.mask.2bit atsei.net2axt &
nohup netToAxt atra.net.netsyn at_vs_ra.prenet arabidopsis_thaliana.mask.2bit rhizophora_apiculata.mask.2bit atra.net2axt &
nohup netToAxt atsa.net.netsyn at_vs_sa.prenet arabidopsis_thaliana.mask.2bit sonneratia_alba.mask.2bit atsa.net2axt &

#b:axtSort
for i in *.net2axt;do axtSort $i $i.sortaxt2;done;

#c:axtToMaf
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=avicennia_marina. atam.net2axt.sortaxt2 arabidopsis_thaliana.chrsize avicennia_marina.chrsize arabidopsis_thaliana.avicennia_marina.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=arabidopsis_thaliana. atat.net2axt.sortaxt2 arabidopsis_thaliana.chrsize arabidopsis_thaliana.chrsize arabidopsis_thaliana.arabidopsis_thaliana.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=eucalyptus_grandis. ateug.net2axt.sortaxt2 arabidopsis_thaliana.chrsize eucalyptus_grandis.chrsize arabidopsis_thaliana.eucalyptus_grandis.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=oryza_sativa. ators.net2axt.sortaxt2 arabidopsis_thaliana.chrsize oryza_sativa.chrsize arabidopsis_thaliana.oryza_sativa.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=populus_trichocarpa. atpot.net2axt.sortaxt2 arabidopsis_thaliana.chrsize populus_trichocarpa.chrsize arabidopsis_thaliana.populus_trichocarpa.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=rhizophora_apiculata. atra.net2axt.sortaxt2 arabidopsis_thaliana.chrsize rhizophora_apiculata.chrsize arabidopsis_thaliana.rhizophora_apiculata.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=sonneratia_alba. atsa.net2axt.sortaxt2 arabidopsis_thaliana.chrsize sonneratia_alba.chrsize arabidopsis_thaliana.sonneratia_alba.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=sesamum_indicum. atsei.net2axt.sortaxt2 arabidopsis_thaliana.chrsize se.indicum.chrsize arabidopsis_thaliana.sesamum_indicum.sing.maf &

#step8:combine multiple alignment
roast + E=arabidopsis_thaliana "(oryza_sativa ((((sonneratia_alba eucalyptus_grandis) arabidopsis_thaliana) (populus_trichocarpa rhizophora_apiculata)) (sesamum_indicum avicennia_marina)))" *.*.sing.maf 7species.final.maf









