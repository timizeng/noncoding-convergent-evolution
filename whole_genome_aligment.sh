#获得基因组序列文件，并修改fa序列名（待补）
#step1: 获得2bit文件，genome.fa ==> genome.2bit, 并提取基因组的 chromesize
#a:写成循环
for i in *.fa;
do
 gzip $i 
 faToTwoBit $i.gz $i.2bit 
 twoBitInfo $i.2bit stdout | sort -k2,2nr > $i.chromsizes
done;
#b:单个运行
gzip ps128.fa 
faToTwoBit ps128.fa.gz ps128.2bit 
twoBitInfo ps128.2bit stdout | sort -k2,2nr > ps128.chromsizes

#c:chromesize文件格式示例：
Chr1	43270923
Chr3	36413819
Chr2	35937250
Chr4	35502694
Chr6	31248787
Chr5	29958434

#step2: reapetmask, 将基因组中的重复序列用 trf 软件进行 mask
for i in *.fa;
do 
 trf $i 2 7 7 80 10 50 500 -f -d -m -h 
done;

#step3: pairwise alignment, 以arabidopsis_thaliana为参考，将其他基因组与之比对，获得pairwise比对结果
for i in *.fa;
do
 lastz-1.04.00 arabidopsis_thaliana.fa[multiple] $i[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=axt > at_vs_$i.axt 
done;

#step4:axt文件转换为psl，并使用axtChain进行chain,例如：
nohup axtToPsl at_vs_at.chain.mask.axt arabidopsis_thaliana.chrsize arabidopsis_thaliana.chrsize at_vs_at.psl &
nohup axtToPsl at_vs_eug.chain.mask.axt arabidopsis_thaliana.chrsize eucalyptus_grandis.chrsize at_vs_eug.psl &
nohup axtToPsl at_vs_ors.chain.mask.axt arabidopsis_thaliana.chrsize oryza_sativa.chrsize at_vs_ors.psl &
nohup axtToPsl at_vs_pot.chain.mask.axt arabidopsis_thaliana.chrsize populus_trichocarpa.chrsize at_vs_pot.psl &
nohup axtToPsl at_vs_ra.chain.mask.axt arabidopsis_thaliana.chrsize rhizophora_apiculata.chrsize at_vs_ra.psl &
nohup axtToPsl at_vs_sa.chain.mask.axt arabidopsis_thaliana.chrsize sonneratia_alba.chrsize at_vs_sa.psl &
nohup axtToPsl at_vs_sei.chain.mask.axt arabidopsis_thaliana.chrsize se.indicum.chrsize at_vs_sei.psl &：
nohup axtToPsl at_vs_am.chain.mask.axt arabidopsis_thaliana.chrsize avicennia_marina.chrsize at_vs_am.psl &

.......将所有文件都转为psl格式......

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
nohup axtToMaf atam.net2axt.sortaxt2 arabidopsis_thaliana.chrsize avicennia_marina.chrsize at_am.maf &
nohup axtToMaf atat.net2axt.sortaxt2 arabidopsis_thaliana.chrsize arabidopsis_thaliana.chrsize at_at.maf &
nohup axtToMaf ateug.net2axt.sortaxt2 arabidopsis_thaliana.chrsize eucalyptus_grandis.chrsize at_eug.maf &
nohup axtToMaf ators.net2axt.sortaxt2 arabidopsis_thaliana.chrsize oryza_sativa.chrsize at_ors.maf &
nohup axtToMaf atpot.net2axt.sortaxt2 arabidopsis_thaliana.chrsize populus_trichocarpa.chrsize at_pot.maf &
nohup axtToMaf atra.net2axt.sortaxt2 arabidopsis_thaliana.chrsize rhizophora_apiculata.chrsize at_ra.maf &
nohup axtToMaf atsa.net2axt.sortaxt2 arabidopsis_thaliana.chrsize sonneratia_alba.chrsize at_sa.maf &
nohup axtToMaf atsei.net2axt.sortaxt2 arabidopsis_thaliana.chrsize se.indicum.chrsize at_sei.maf &


#axtToMaf可以改 maf中的序列标签，-tPrefix 和 -qPrefix，最后的roast运行,要求 maf文件中的比对标签必须带上物种名,且与输入物种树中的物种名完全一致
axtToMaf -tPrefix=$tn -qPrefix=$qn $output_dir/6.net_to_axt/all_sort.axt $output_dir/$tn.sizes $output_dir/$qn.sizes $output_dir/7.maf/all.maf

#dtplot
for i in *.sing.maf;do last-dotplot $i $i.png;done;

#step8:combine multiple alignment
#maf内的比对标签，输入文件名，需要与输入物种树名称完全一致
roast + E=arabidopsis_thaliana "(oryza_sativa ((((sonneratia_alba eucalyptus_grandis) arabidopsis_thaliana) (populus_trichocarpa rhizophora_apiculata)) (sesamum_indicum avicennia_marina)))" *.*.sing.maf 7species.final.maf

#phyloFit 获得nonconserved.mod
phyloFit --tree "(oryza_sativa,((((sonneratia_alba,eucalyptus_grandis),arabidopsis_thaliana),(populus_trichocarpa,rhizophora_apiculata)),(sesamum_indicum,avicennia_marina)))" --msa-format MAF --out-root nonconserved-4d 7species.final.maf


#使用 mafSplit 分解比对文件
mafSplit -byTarget arabidopsis_thaliana -useFullSequenceName outdir/ 7species.final.maf