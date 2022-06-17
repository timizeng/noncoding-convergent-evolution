#获得基因组序列文件，并修改fa序列名（待补），注意不要包含 | 等多余字符，可以有空格，要符合chrsize格式
#step0:去除重复序列，reaptemask
nohup trf /public1/users/zengws/workspace/convergent-evo/cnnes/whole-genome-aliment/geomedb/avicennia_marina/AM_final.fa 2 7 7 80 10 50 500 -f -d -m -h &
nohup trf /public1/users/zengws/workspace/convergent-evo/cnnes/whole-genome-aliment/geomedb/oryza_sativa/osativa_v7.fa 2 7 7 80 10 50 500 -f -d -m -h &
nohup trf /public1/users/zengws/workspace/convergent-evo/cnnes/whole-genome-aliment/geomedb/e_grandis/egrandis_v2.fa 2 7 7 80 10 50 500 -f -d -m -h &
nohup trf /public1/users/zengws/workspace/convergent-evo/cnnes/whole-genome-aliment/geomedb/rhizophora_apiculata/RA_final.fa 2 7 7 80 10 50 500 -f -d -m -h &
nohup trf /public1/users/zengws/workspace/convergent-evo/cnnes/whole-genome-aliment/geomedb/populus_trichocarpa/ptrichocarpa_v4.fa 2 7 7 80 10 50 500 -f -d -m -h &
nohup trf /public1/users/zengws/workspace/convergent-evo/cnnes/whole-genome-aliment/geomedb/sonneratia_alba/SA_final.fa 2 7 7 80 10 50 500 -f -d -m -h &
nohup trf /public1/users/zengws/workspace/convergent-evo/cnnes/whole-genome-aliment/geomedb/se_indicum/se.indicum.dna.fa 2 7 7 80 10 50 500 -f -d -m -h &

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
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=avicennia_marina. atam.net2axt.sortaxt2 arabidopsis_thaliana.chrsize avicennia_marina.chrsize arabidopsis_thaliana.avicennia_marina.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=arabidopsis_thaliana. atat.net2axt.sortaxt2 arabidopsis_thaliana.chrsize arabidopsis_thaliana.chrsize arabidopsis_thaliana.arabidopsis_thaliana.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=eucalyptus_grandis. ateug.net2axt.sortaxt2 arabidopsis_thaliana.chrsize eucalyptus_grandis.chrsize arabidopsis_thaliana.eucalyptus_grandis.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=oryza_sativa. ators.net2axt.sortaxt2 arabidopsis_thaliana.chrsize oryza_sativa.chrsize arabidopsis_thaliana.oryza_sativa.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=populus_trichocarpa. atpot.net2axt.sortaxt2 arabidopsis_thaliana.chrsize populus_trichocarpa.chrsize arabidopsis_thaliana.populus_trichocarpa.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=rhizophora_apiculata. atra.net2axt.sortaxt2 arabidopsis_thaliana.chrsize rhizophora_apiculata.chrsize arabidopsis_thaliana.rhizophora_apiculata.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=sonneratia_alba. atsa.net2axt.sortaxt2 arabidopsis_thaliana.chrsize sonneratia_alba.chrsize arabidopsis_thaliana.sonneratia_alba.sing.maf &
nohup axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=sesamum_indicum. atsei.net2axt.sortaxt2 arabidopsis_thaliana.chrsize se.indicum.chrsize arabidopsis_thaliana.sesamum_indicum.sing.maf &


#axtToMaf可以改 maf中的序列标签，-tPrefix 和 -qPrefix，最后的roast运行,要求 maf文件中的比对标签必须带上物种名,且与输入物种树中的物种名完全一致
axtToMaf -tPrefix=$tn -qPrefix=$qn $output_dir/6.net_to_axt/all_sort.axt $output_dir/$tn.sizes $output_dir/$qn.sizes $output_dir/7.maf/all.maf

#dtplot
for i in *.sing.maf;do last-dotplot $i $i.png;done;

#step8:combine multiple alignment
#maf内的比对标签，输入文件名，需要与输入物种树名称完全一致
roast + E=arabidopsis_thaliana "(oryza_sativa ((((sonneratia_alba eucalyptus_grandis) arabidopsis_thaliana) (populus_trichocarpa rhizophora_apiculata)) (sesamum_indicum avicennia_marina)))" *.*.sing.maf 7species.final.maf

#使用 mafSplit 分解比对文件,以arabidopsis_thaliana染色体为单位分解
mafSplit -byTarget arabidopsis_thaliana -useFullSequenceName mafplit/ 7species.final.maf

#phyloFit提取4d位点构建,获得nonconserved.mod
#step1：从gtf文件中提取各自染色体的注释
less -S  Araport11.gtf | grep 'ChrM' | grep 'CDS'> ChrM.gtf
#step2:替换gtf中的染色体名使之与全基因组比对中的maf文件中的名称一致
sed -i 's/ChrM/arabidopsis_thaliana.ChrM/' ChrM.gtf
#step3：使用msa_view提取4dcodon位点
msa_view --in-format MAF chrM.maf --4d --features ChrM.gtf --out-format SS > ChrM.4dcodon.ss
#step4:使用msa_view 提取4d site
msa_view ChrM.4dcodon.ss --in-format SS --out-format SS --tuple-size 1 > ChrM.4dsites.ss

#写成循环
for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; do cat Araport11.gtf | grep $i | grep 'CDS' > $i.gtf; done;
sed -i 's/Chr1/arabidopsis_thaliana.Chr1/' Chr1.gtf
.....全部替换，不知为何 sed 无法允许循环.....
for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; do msa_view $i.maf --4d --features $i.gtf > $i.4d-codons.ss ; done;
for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; do msa_view $i.4d-codons.ss --in-format SS --out-format SS --tuple-size 1 > $i.4d-sites.ss; done;

#step5:使用 msa_view整合每个染色体的4dsites.
msa_view --unordered-ss --out-format SS --aggregate arabidopsis_thaliana,populus_trichocarpa,sonneratia_alba,oryza_sativa,eucalyptus_grandis,sesamum_indicum,avicennia_marina,rhizophora_apiculata *.4d-sites.ss > all.4d-site.ss
#phyloFit生成 nonconserved-4d.mod
phyloFit --tree "(oryza_sativa,((((sonneratia_alba,eucalyptus_grandis),arabidopsis_thaliana),(populus_trichocarpa,rhizophora_apiculata)),(sesamum_indicum,avicennia_marina)))" --msa-format SS --out-root nonconserved-4d all.4d-site.ss


#使用msa_split分解单条染色体比对的maf
msa_split Chr1.maf --in-format MAF --windows 100000,0 --out-root split/Chr1 --out-format SS --min-informative 1000 --between-blocks 5000
msa_split Chr2.maf --in-format MAF --windows 100000,0 --out-root split/Chr2 --out-format SS --min-informative 1000 --between-blocks 5000
msa_split Chr3.maf --in-format MAF --windows 100000,0 --out-root split/Chr3 --out-format SS --min-informative 1000 --between-blocks 5000
msa_split Chr4.maf --in-format MAF --windows 100000,0 --out-root split/Chr4 --out-format SS --min-informative 1000 --between-blocks 5000
msa_split Chr5.maf --in-format MAF --windows 100000,0 --out-root split/Chr5 --out-format SS --min-informative 1000 --between-blocks 5000
msa_split ChrC.maf --in-format MAF --windows 100000,0 --out-root split/ChrC --out-format SS --min-informative 1000 --between-blocks 5000
msa_split ChrM.maf --in-format MAF --windows 100000,0 --out-root split/ChrM --out-format SS --min-informative 1000 --between-blocks 5000

for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; do msa_split $i.maf --in-format MAF --windows 100000,0 --out-root windowsplit/$i --out-format SS --min-informative 1000 --between-blocks 5000; done
#run phastcons, 使用生成的.ss文件，运行phastcons生成 conserved.mod 和 nonconserverd.mod
for i in *.ss;
do
        phastCons --estimate-rho $i --no-post-probs $i /data/users/zengws/cnees/test/rpmask-genome/mafsplit_file/nonconserved-4d.mod
done;

输出一系列的cons.mod和noncons.mod文件，cons.mod用于下一步生成Average.chr.cons.mod
#一下3个ss在计算mod时卡死，删除后计算
chr4.18499916-18583278.ss
chr5.26900001-26965661.ss
chrc.99627-154133.ss
#run phyloBoot, 对每一条分解后的染色体运行，生成该条染色体的保守模型的平均模型 Average.chr.cons.mod
phyloBoot --read-mods *.cons.mod --output-average Average.chr.cons.mod、

其中输入文件为window.cons.list, 该list为每个mod文件名，由 , 隔开；
输出文件为Average.chr.cons.mod

#run phastCons, 预测cnees，并输出预测cnees位置的bed文件，以参考基因组arabidopsis_thaliana的位置输出
nohup phastCons --most-conserved phastCons.Chr1.bed --score Chr1.maf Average.Chr1.cons.mod,nonconserved-4d.mod > phastCons.Chr1.wig &
nohup phastCons --most-conserved phastCons.Chr2.bed --score Chr2.maf Average.Chr2.cons.mod,nonconserved-4d.mod > phastCons.Chr2.wig &
nohup phastCons --most-conserved phastCons.Chr3.bed --score Chr3.maf Average.Chr3.cons.mod,nonconserved-4d.mod > phastCons.Chr3.wig &
nohup phastCons --most-conserved phastCons.Chr4.bed --score Chr4.maf Average.Chr4.cons.mod,nonconserved-4d.mod > phastCons.Chr4.wig &
nohup phastCons --most-conserved phastCons.Chr5.bed --score Chr5.maf Average.Chr5.cons.mod,nonconserved-4d.mod > phastCons.Chr5.wig &
nohup phastCons --most-conserved phastCons.ChrM.bed --score ChrM.maf Average.ChrM.cons.mod,nonconserved-4d.mod > phastCons.ChrM.wig &
nohup phastCons --most-conserved phastCons.ChrC.bed --score ChrC.maf Average.ChrC.cons.mod,nonconserved-4d.mod > phastCons.ChrC.wig &
