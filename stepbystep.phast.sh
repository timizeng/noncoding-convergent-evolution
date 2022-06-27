#step1 长度筛选，只保留基因组中大于10000bp的序列
for i in *.fa;
do 
    base=$(basename ${i} .fa)
    echo "nohup seqkit seq -m 10000 --remove-gaps ${i} > ${base}.lenfil.fa &" >>length_fil.sh
done

chmod 755 length_fil.sh
./length_fil.sh


#step2 trf
for i in *.lenfil.fa;
do 
    base=$(basename ${i} .lenfil.fa) 
    echo "nohup trf ${i} 2 7 7 80 10 50 500 -f -d -m -h &" >> trf.sh
done

chmod 755 trf.sh
./trf.sh

for i in *.mask; do base=$(basename ${i} .lenfil.fa.2.7.7.80.10.50.500.mask); mv ${i} ${base}.mask.fa;done;
rm *.dat

#step3 fatotwobit and chromsize
for i in *.fa; 
do 
    base=$(basename ${i} .mask.fa)
    faToTwoBit ${i} ${base}.2bit
    twoBitInfo ${base}.2bit stdout | sort -k2,2nr > ${base}.chromsizes
done

#lastz
for i in *.mask.fa;
do 
    base=$(basename ${i} .mask.fa) 
    echo "nohup lastz-1.04.00 arabidopsis_thaliana.mask.fa[multiple] \
    ${i}[multiple] --strand=both --notransition --ambiguous=iupac --step=20 \
    --nogapped --format=axt > ${base}.axt &" >> lastz.sh
done
chmod 755 lastz.sh
./lastz.sh

#axt to psl
for i in *.axt;
do 
    base=$(basename ${i} .axt)
    axtToPsl ${i} ../chromsize/arabidopsis_thaliana.chromsizes ../chromsize/${base}.chromsizes ${base}.psl
done

#axtChain
for i in *.psl; 
do 
    base=$(basename ${i} .psl)
    axtChain -linearGap=loose -psl ${i} \
    -faT ../trfgenome/arabidopsis_thaliana.mask.fa \
    -faQ ../trfgenome/${base}.mask.fa ${base}.chain
done

#Netting,chainNet
##a:chainMergeSort
for i in *.chain;
do
    base=$(basename ${i} .chain)
	chainMergeSort ${i} > ${base}.chain.sort
done

#b:chainPreNet
for i in *.chain.sort;
do
    base=$(basename ${i} .chain.sort)
	chainPreNet ${i} ../chromsize/arabidopsis_thaliana.chromsizes \
    ../chromsize/${base}.chromsizes ${base}.prenet
done

#c:chainNet
for i in *.prenet;
do
    base=$(basename ${i} .prenet)
    chainNet ${i} ../chromsize/arabidopsis_thaliana.chromsizes \
    ../chromsize/${base}.chromsizes ${base}.at.net ${base}.net
done

#d:netSyntenic
#d:netSyntenic
for i in *.at.net
do 
    base=$(basename ${i} .at.net)
    netSyntenic ${i} ${base}.netsyn
done

#maffting
#a:netToAxt
for i in *.netsyn
do 
    base=$(basename ${i} .netsyn)
    netToAxt ${i} ${base}.prenet ../twobitgenome/arabidopsis_thaliana.2bit \
    ../twobitgenome/${base}.2bit ${base}.net2axt
done

#b:axtSort
for i in *.net2axt;
do 
    base=$(basename ${i} .net2axt)
    axtSort ${i} ${base}.sort.axt
done

#c:axtToMaf
for i in *.sort.axt;
do 
    base=$(basename ${i} .sort.axt)
    axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=${base}. \
    ${i} ../chromsize/arabidopsis_thaliana.chromsizes \
    ../chromsize/${base}.chromsizes ${base}.arabidopsis_thaliana.sing.maf 
done

#重命名一下，必须要把arabidopsis_thaliana在前
for i in *.sing.maf; do base=$(basename ${i} .arabidopsis_thaliana.sing.maf); mv ${i} arabidopsis_thaliana.${base}.sing.maf; done;

#step8:combine multiple alignment
roast + E=arabidopsis_thaliana "(((arabidopsis_thaliana theobroma_cacao) (citrus_sinensis (swietenia_macrophylla (xylocarpus_granatum (xylocarpus_rumphii xylocarpus_moluccensis))))) (populus_trichocarpa (carallia_pectinifolia ((bruguiera_gymnorhiza bruguiera_sexangula) ((ceriops_tagal (kandelia_obovata kandelia_candel)) (rhizophora_mangle (rhizophora_apiculata (rhizophora_mucronata rhizophora_stylosa))))))))" *.*.sing.maf 18species.final.maf

.................................以下开始phast流程.............................
####phast
#使用 mafSplit 分解比对文件,以arabidopsis_thaliana染色体为单位分解
mafSplit -byTarget arabidopsis_thaliana -useFullSequenceName mafplit/ 18species.final.maf


...........................提取4d位点估计单条染色体non.conserved.mod.....................
#提取4d位点估计 non.conserved.mod
for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; 
do 
    msa_view ${i}.maf --in-format MAF --4d --features ${i}.gtf > 4dsite/${i}.4d-codons.ss 
done

for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; 
do 
    msa_view ${i}.4d-codons.ss --in-format SS --out-format SS --tuple-size 1 > ${i}.4d-sites.ss
done

#合并所有染色体4d位点
msa_view --unordered-ss --out-format SS --aggregate arabidopsis_thaliana,theobroma_cacao,citrus_sinensis,swietenia_macrophylla,xylocarpus_granatum,xylocarpus_rumphii,xylocarpus_moluccensis,populus_trichocarpa,carallia_pectinifolia,bruguiera_gymnorhiza,bruguiera_sexangula,ceriops_tagal,kandelia_obovata,kandelia_candel,rhizophora_mangle,rhizophora_apiculata,rhizophora_mucronata,rhizophora_stylosa *sites.ss > all.4d-sites.ss

#估计总体的non.cons.mod
phyloFit --tree "(((arabidopsis_thaliana,theobroma_cacao),(citrus_sinensis,(swietenia_macrophylla,(xylocarpus_granatum,(xylocarpus_rumphii,xylocarpus_moluccensis))))),(populus_trichocarpa,(carallia_pectinifolia,((bruguiera_gymnorhiza,bruguiera_sexangula),((ceriops_tagal,(kandelia_obovata,kandelia_candel)),(rhizophora_mangle,(rhizophora_apiculata,(rhizophora_mucronata,rhizophora_stylosa))))))))" --msa-format SS --out-root nonconserved-4d all.4d-sites.ss

#计算单条染色体的non.cons.mod
for i in *.4d-sites.ss;
do 
    base=$(basename ${i} .4d-sites.ss)
    phyloFit --tree "(((arabidopsis_thaliana,theobroma_cacao),(citrus_sinensis,(swietenia_macrophylla,(xylocarpus_granatum,(xylocarpus_rumphii,xylocarpus_moluccensis))))),(populus_trichocarpa,(carallia_pectinifolia,((bruguiera_gymnorhiza,bruguiera_sexangula),((ceriops_tagal,(kandelia_obovata,kandelia_candel)),(rhizophora_mangle,(rhizophora_apiculata,(rhizophora_mucronata,rhizophora_stylosa))))))))" --msa-format SS --out-root nonconserved-4d.${base} ${i}
done


.......................通过estimate-rho分别计算单条染色体的con.elements....................
##通过estimate-rho和每条染色体4d位点估计出的noncons.mod来估算每条染色体的cons.mod
for i in *.maf;
do 
    base=$(basename ${i} .maf)
    echo "nohup phastCons --target-coverage 0.25 --expected-length 12 \
    --rho 0.4 --estimate-rho ${base}.rho --msa-format MAF \
    ${base}.maf 4dsite/nonconserved-4d.${base}.mod \
    --no-post-probs &" >> estimate-rho.sh
done
chmod 755 estimate-rho.sh
./estimate-rho.sh

#计算根据每条染色体的cons.mod和non.cons.mod计算conseverd.element
for i in *.maf; 
do      
    base=$(basename ${i} .maf)   
    echo "nohup phastCons --target-coverage 0.25 --expected-length 12 \
    --rho 0.4 --msa-format MAF ${base}.maf --most-conserved phastCons.${base}.bed \
    ${base}.rho.cons.mod,${base}.rho.noncons.mod > phastCons.${base}.wig &" >> con.element.sh 
done
chmod 755 con.element.sh
./con.element.sh

......该方法出现警告.....如下：
Finding MLE for (rho)...
Computing emission probs (state 1, cat 1, mod 1)...
WARNING: likelihood decreased during EM: it 8 total_logl=-580783.8808, prev_total_logl=-580783.8801
(rho = 0.199658)


.................通过滑窗和 --estimate-trees 计算滑窗的cons.mod 和 non.cons.mod 最后对单个染色体聚合mods......................
##通过滑窗计算 每条 染色体的 cons.mod
#根据 http://compgen.cshl.edu/phast/phastCons-HOWTO.html 中的4.1部分教程重新分析
1.Split the alignments into small fragments,当滑窗窗口变小时，
--min-informative 1000 --between-blocks 5000 需要适当变小，
不然会丢弃掉很多block

mkdir -p chunks            # put fragments here
        for file in Chr*.maf ; do
            root=$(basename ${file} .maf)
            msa_split ${file} --in-format MAF --refseq ${root}.fa \
                --windows 1000000,0 --out-root chunks/${root} --out-format SS \
                --min-informative 1000 --between-blocks 5000 
        done

2.Estimate parameters for each fragment.

mkdir -p trees     # put estimated tree models here
rm -f trees/*      # in case old versions left over

#以4d位点得到的对应染色体non.cons.mod作为输入，每个染色体逐条运行
for file in chunks/Chr1.*.ss ; 
do 
    root=$(basename ${file} .ss)  
    phastCons --target-coverage 0.125 --expected-length 20 \
    --gc 0.4 --estimate-trees trees/${root} ${file} \
    nonconserved-4d.Chr1.mod --no-post-probs
done

3.Combine the separately estimated parameters by phyloBoot

for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC;
do
    ls ../trees/${i}.*.cons.mod > ${i}.all.cons.txt
    ls ../trees/${i}.*.noncons.mod > ${i}.all.noncons.txt 
done

for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; #此处要逐个运行，循环有问题
do
    phyloBoot --read-mods '*${i}.all.cons.txt' --output-average ${i}.ave.cons.mod
    phyloBoot --read-mods '*${i}.all.noncons.txt' --output-average ${i}.ave.noncons.mod
done


4.Predict conserved elements and conservation scores globally using the combined estimates
#对每条染色体单独运行
for file in chunks/Chr1.*.ss ;  
do 
    root=$(basename ${file} .ss)
    phastCons --target-coverage 0.125 --expected-length 20 \
    --most-conserved elements/${root}.bed --score ${file} \
    averagemod/Chr1.ave.cons.mod,averagemod/Chr1.ave.noncons.mod > scores/${root}.wig
done


....................评价non.cons.mod.........
#分割chunk, --min-informative 1000 --between-blocks 5000 要根据 window大小适当缩放
for file in Chr*.maf ;
do             
    root=$(basename ${file} .maf)             
    msa_split ${file} --in-format MAF --refseq ${root}.fa --windows 1000000,0\
     --out-root chunks/${root} --out-format SS --min-informative 1000 \
     --between-blocks 5000
done

#将chunk分到多个文件 Chr1.1 Chr1.2 ...... Chr5.5 分成子块运行，每条染色体分5个文件
#创建多个文件，将对应的chunk转移到对应文件夹运行估算mod过程，
for i in 1 2 3 4 5; do for j in 1 2 3 4 5 ; do mkdir -p Chr${i}.${j};done;done;
#nohup挂起滑窗mod估计过程
nohup sh -c 'for file in Chr5.5/Chr5.*.ss;  do root=$(basename ${file} .ss)  ; phastCons --estimate-trees trees/${root} ${file} all.nonconserved-4d.mod --no-post-probs; done; ' &

#从mod中提取数字
grep "TREE" ChrM.300001-366651.cons.mod |sed -e 's/[A-Za-z : _ () ; ,]/ /g' | tr -s ' ' | sed 's/ /\n/g'

sed -e 's/[A-Za-z : _ () ; ,]/ /g' 将任何大小写字母 : _ () ; , 都替换成空格
tr -s ' ' 将多个空格压缩成单个空格
sed 's/ /\n/g' 空格之间换行

写成for循环
for i in Chr*.*.mod; do grep "TREE" ${i} |sed -e 's/[A-Za-z : _ () ; ,]/ /g' | tr -s ' ' | sed 's/ /\n/g' > ${i}.txt ; done;

再筛选计算
for i in Chr*.*.noncons.mod.txt; do cat ${i} | awk '$1>=0 && $1<=1 {print $1}' | awk '{sum+=$1}END{print sum}' >> noncons.final.txt;done;


一步到位
for i in Chr*.*.noncons.mod;do grep "TREE" ${i} | sed -e 's/[A-Za-z : _ () ; ,]/ /g' | tr -s ' ' | sed 's/ /\n/g' |  awk '$1>=0 && $1<=1 {print $1}' | awk '{sum+=$1}END{print sum}' >> noncons.txt;done;


all.nonconserved-4d.mod 支长和为 3.87056


............是用repeat估计non.cons.mod.................

phyloFit --tree "(((arabidopsis_thaliana,theobroma_cacao),(citrus_sinensis,(swietenia_macrophylla,(xylocarpus_granatum,(xylocarpus_rumphii,xylocarpus_moluccensis))))),(populus_trichocarpa,(carallia_pectinifolia,((bruguiera_gymnorhiza,bruguiera_sexangula),((ceriops_tagal,(kandelia_obovata,kandelia_candel)),(rhizophora_mangle,(rhizophora_apiculata,(rhizophora_mucronata,rhizophora_stylosa))))))))" 
--features Chr1.gtf --do-cats repeat --msa-format MAF --out-root Chr1.repeat.non.cons Chr1.maf


nohup sh -c 'for file in Chr1.1/Chr1.*.ss;  do root=$(basename ${file} .ss)  ; phastCons --estimate-trees LINEtrees/${root} ${file} ../Average.LINE.non.cons.mod --no-post-probs; done; ' &



........在不同物种之间转换坐标.....liftover可行.....
liftOver conserved.bed at_vs_ors.chain orsconserved.bed unmapp.bed
#liftover中有一个参数 -minMatch=0.95 (默认值)， 可以调整，使得比对更严格或者更宽松 
#                    -multiple 设置允许输出多重比对，如果不设置，则多比对的位置都不输出（建议设置）


.......使用mafFind也可以做.........例如有一个pair-wirse对齐块：

mafFind arabidopsis_thaliana.bruguiera_gymnorhiza.sing.maf 41693 41708 arabidopsis_thaliana.ChrM slice > out.file


再从中提取坐标，若目标物种对应到的是负链 (-)，则转换坐标的方式为：
a score=51656.000000
s arabidopsis_thaliana.ChrM                     41693 15 +   366924 ACCACTCGCTCGCCA
s bruguiera_gymnorhiza.unplaced_scaffold1080     8235 15 -    10476 ACCACTCGCTCGCCG

    染色体位置             开始位置         结束位置     名字（随便）  分数（随便）   链（重要）
unplaced_scaffold1080   10476-8235-15    10476-8235        name1       0           -

再用bedtools getfasta 指定 -s 提取目标序列

若目标物种对应到的是正链 (+)，则转换坐标的方式为：
s arabidopsis_thaliana.ChrC                112085 250 + 154478 CTACACGATTAGTTACAAAAGCTTTTTGACAAGCATTCGCTGCAATAGGTCGTGTGAACCAAAAACCTATTAATAAATACGAACACATTCCAACTAATTCCCAAAAAAAATAAACTTGGATCAAATTAGAACTAGTAACTAATCCTAACATTGAAGTATT--AAAAAAACCCATATAAGCAAAAAATCTCAGGTATCCTTGATCATGCGACATATAATTGTCACTATAAATCAGAACCAAAATCCCAACAGT
s bruguiera_gymnorhiza.unplaced_scaffold63  26350 252 +  27206 CTACACGATTAGTTATAAATGCCTTTTGACAAGCATTCGAGGCAATAGGTTGTGTGAACCAAAACCCTATTAATAAATAAGAACACATTCCAACCAATTCCCAAAAAATATAAATTTGTATCAAATTAGAACTAGTAACTAATCCCACCATTGAAGTATTGAAAAAAAAATCATATAAGCAAAAAATCTCAAATAACTTTGATCATGAGCCATATAATTGTCACTATAAAGAAGAACCAAAGTTCCAACTGT

染色体位置             开始位置    结束位置     名字（随便） 分数（随便）  链（重要）
unplaced_scaffold63    26350     26350+252      name1       0           +
