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

phastCons --most-conserved phastCons.chr.bed --score Chr.maf Average.chr.cons.mod,nonconserved-4d.mod > phastCons.chr.wig
输入文件为mafSplit步根据染色体分解的Chr.maf，phyloBoot步的Average.chr.cons.mod，以及phyloFit步的nonconserved-4d.mod
输出文件为phastCons.chr.bed（预测得到的cnees位置）和phastCons.chr.wig（）
bed文件：
染色体  strat    end      序号   得分    链
ChrC    410     496     ChrC.1  38      +
ChrC    532     1405    ChrC.2  372     +
ChrC    1516    1564    ChrC.3  26      +
ChrC    1716    1874    ChrC.4  54      +
ChrC    2227    2301    ChrC.5  31      +
ChrC    2645    2768    ChrC.6  17      +
ChrC    4192    4254    ChrC.7  33      +
ChrC    4277    4347    ChrC.8  38      +
ChrC    6136    6205    ChrC.9  34      +

#将每个染色体预测得到的保守元件bed文件合成一个
cat *.bed | sort -k1,1 -k2,2n > conserved.bed
awk '{print NR}' conserved.bed | tail -n1 #3059
#统计cnees长度分布并作图 R
cat conserved.bed | awk '{print $3-$2}' > cneelength.txt

library(tidyverse)
length <- read.table(file = 'cneelength.txt', header = F)
length[length>1000]=1001 #将大于1000的值变为1001，便于统计

breaks <- seq(0,1050,50)
data <- as.data.frame(table(cut(length$V1,breaks))) #分区间统计
data <- edit(data) #把最后一个区间标签手动换成 >1000 

ggplot(data = data,aes(x=data$Var1,y=data$Freq))+
  geom_col(fill="#69b3a2",color="#e9ecef")+
  labs(x='length distribution', y='conserved elements count')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size=8))+
  geom_text(aes(label=data$Freq),size=3, vjust=-0.5)

#统计转录起始位点TSS区域富集情况

##根据bed文件，删除长度小于等于50bp的element
cat conserved.bed | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$3-$2}' > cons.bed
cat cons.bed | awk -v OFS='\t' '{if($7>50) print $1,$2,$3,$4,$5,$6 }' > cons.filter50.bed

#根据arabidopsis_thaliana的gff文件，与上一步得到的bed文件，删除与exon,CDS等有overlapping的区域
#a:检查拟南芥gff文件，该文件只包含exon和cds,因此不需要处理直接
#b:使用bedtools分析overlap
bedtools intersect -wao -a cons.filter50.bed -b refGene.exons.gtf > overlap.bed

#c：查看overlap.bed文件，删除和exon和CDS有overlap的行
Chr1    111889  111961  Chr1.2  72      +       .       .       .       -1      -1      .       .       .       .       0
Chr1    309274  309351  Chr1.8  120     +       .       .       .       -1      -1      .       .       .       .       0
Chr1    515493  515567  Chr1.9  126     +       .       .       .       -1      -1      .       .       .       .       0
Chr1    552639  552711  Chr1.11 123     +       Chr1    phytozomev12    CDS     552687  552761  .       -       0       transcript_id "AT1G02590.1.Ara
Chr1    608542  608610  Chr1.13 32      +       Chr1    phytozomev12    CDS     608398  608711  .       -       2       transcript_id "AT1G02780.1.Ara
Chr1    741681  741754  Chr1.14 28      +       Chr1    phytozomev12    exon    741684  741795  .       +       .       transcript_id "AT1G03090.2.Ara
Chr1    741681  741754  Chr1.14 28      +       Chr1    phytozomev12    CDS     741684  741795  .       +       0       transcript_id "AT1G03090.2.Ara
Chr1    741681  741754  Chr1.14 28      +       Chr1    phytozomev12    exon    741684  741795  .       +       .       transcript_id "AT1G03090.1.Ara
Chr1    741681  741754  Chr1.14 28      +       Chr1    phytozomev12    CDS     741684  741795  .       +       0       transcript_id "AT1G03090.1.Ara
#没有overlap的第10，11列为 -1 ，输出包含 -1 列 的行，只需要前6列就行
cat overlap.bed | awk -v OFS='\t' '{if($10==-1) print $1,$2,$3,$4,$5,$6}' > nonexon.element.bed
awk '{print NR}' nonexon.element.bed | tail -n1 #1297个

#再次统计长度分布以及总数



#例子，对于A文件中染色体位置，如果和B文件中染色体位置有overlap,
#则输出在A文件中染色体位置和在B文件中染色体位置，以及overlap的长度；
#如果和B文件中染色体位置都没有overlap,则用'. -1-1'补齐文件
$ cat A.bed 
chr1 10 20 
chr1 30 40 
$ cat B.bed
chr1 15 20
chr1 18 25
$ bedtools intersect -a A.bed -b B.bed -wao
chr1 10 20 chr1 15 20 5
chr1 10 20 chr1 18 25 2
chr1 30 40 . -1 -1


#使用Lifeover转换坐标，本质就是Blastn
liftOver conserved.bed at_vs_ors.chain orsconserved.bed unmapp.bed
#liftover中有一个参数 -minMatch=0.95 (默认值)， 可以调整，使得比对更严格或者更宽松 
#                    -mulitple 设置允许输出多重比对，如果不设置，则多比对的位置都不输出（建议设置）
#bedtools提取对应序列
bedtools getfasta -fi oryza_sativa.mask.fa -bed orsconserved.bed -s > orycnees.fa



ALPHABET: A C G T
ORDER: 0
SUBST_MOD: REV
TRAINING_LNL: -44411214.590680
BACKGROUND: 0.277898 0.222206 0.221925 0.277971
RATE_MAT:
  -0.959657    0.197598    0.491817    0.270242
   0.247122   -1.045227    0.182576    0.615528
   0.615862    0.182808   -1.050870    0.252200
   0.270171    0.492045    0.201350   -0.963566
TREE: (oryza_sativa:0.109119,((((sonneratia_alba:0.141679,eucalyptus_grandis:0.153766):0.0424577,arabidopsis_thaliana:0.176973):0.0191425,(populus_trichocarpa:0.133233,rhizophora_apiculata:0.132441):0.0402056):0.0244227,(sesamum_indicum:0.120677,avicennia_marina:0.130972):0.0516579):0.109119);


ALPHABET: A C G T
ORDER: 0
SUBST_MOD: REV
TRAINING_LNL: -4180265.866951
BACKGROUND: 0.258506 0.199520 0.171398 0.370576
RATE_MAT:
  -1.254595    0.284576    0.072066    0.897953
   0.368707   -0.713224    0.194462    0.150054
   0.108691    0.226369   -0.959001    0.623942
   0.626391    0.080790    0.288584   -0.995765
TREE: (oryza_sativa:529.654,((((sonneratia_alba:1854.14,eucalyptus_grandis:1756.48):2024.8,arabidopsis_thaliana:4295.12):1634.22,(populus_trichocarpa:1507.75,rhizophora_apiculata:1307.06):1260.47):952.635,(sesamum_indicum:976.937,avicennia_marina:1028.46):934.912):529.654);

#按条拆分染色体
 awk '/^>/{f=++d".fa"} {print > f}' arabidopsis_thaliana.mask.fa


 #不使用4d位点直接生成 neutral model
 phyloFit --tree "(oryza_sativa,((((sonneratia_alba,eucalyptus_grandis),arabidopsis_thaliana),(populus_trichocarpa,rhizophora_apiculata)),(sesamum_indicum,avicennia_marina)))"  --msa-format MAF --out-root neutral all.final.maf

#phylop
 phyloP --wig-scores --mode CONACC --method LRT neutral.mod chr.maf > phyloP.chr.wig