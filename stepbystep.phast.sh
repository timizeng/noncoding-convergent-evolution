#一键生成脚本
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
    ${i}[multiple] --notransition --ambiguous=iupac --step=20 \
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

####phast
#使用 mafSplit 分解比对文件,以arabidopsis_thaliana染色体为单位分解
mafSplit -byTarget arabidopsis_thaliana -useFullSequenceName mafplit/ 18species.final.maf

#提取4d位点估计 non.conserved.mod
for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; 
do 
    msa_view ${i}.maf --in-format MAF --4d --features ${i}.gtf > 4dsite/${i}.4d-codons.ss 
done

for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; 
do 
    msa_view ${i}.4d-codons.ss --in-format SS --out-format SS --tuple-size 1 > ${i}.4d-sites.ss
done

#
msa_view --unordered-ss --out-format SS --aggregate arabidopsis_thaliana,theobroma_cacao,citrus_sinensis,swietenia_macrophylla,xylocarpus_granatum,xylocarpus_rumphii,xylocarpus_moluccensis,populus_trichocarpa,carallia_pectinifolia,bruguiera_gymnorhiza,bruguiera_sexangula,ceriops_tagal,kandelia_obovata,kandelia_candel,rhizophora_mangle,rhizophora_apiculata,rhizophora_mucronata,rhizophora_stylosa *sites.ss > all.4d-sites.ss

#总体的non.cons.mod
phyloFit --tree "(((arabidopsis_thaliana,theobroma_cacao),(citrus_sinensis,(swietenia_macrophylla,(xylocarpus_granatum,(xylocarpus_rumphii,xylocarpus_moluccensis))))),(populus_trichocarpa,(carallia_pectinifolia,((bruguiera_gymnorhiza,bruguiera_sexangula),((ceriops_tagal,(kandelia_obovata,kandelia_candel)),(rhizophora_mangle,(rhizophora_apiculata,(rhizophora_mucronata,rhizophora_stylosa))))))))" --msa-format SS --out-root nonconserved-4d all.4d-sites.ss

#计算单条染色体的non.cons.mod
for i in *.4d-sites.ss;
do 
    base=$(basename ${i} .4d-sites.ss)
    phyloFit --tree "(((arabidopsis_thaliana,theobroma_cacao),(citrus_sinensis,(swietenia_macrophylla,(xylocarpus_granatum,(xylocarpus_rumphii,xylocarpus_moluccensis))))),(populus_trichocarpa,(carallia_pectinifolia,((bruguiera_gymnorhiza,bruguiera_sexangula),((ceriops_tagal,(kandelia_obovata,kandelia_candel)),(rhizophora_mangle,(rhizophora_apiculata,(rhizophora_mucronata,rhizophora_stylosa))))))))" --msa-format SS --out-root nonconserved-4d.${base} ${i}
done


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


##通过滑窗计算 每条 染色体的 cons.mod
#根据 http://compgen.cshl.edu/phast/phastCons-HOWTO.html 中的4.1部分教程重新分析
1.Split the alignments into small fragments

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

#每个染色体逐条运行
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



