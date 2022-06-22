

#使用 mafSplit 分解比对文件,以arabidopsis_thaliana染色体为单位分解
mafSplit -byTarget arabidopsis_thaliana -useFullSequenceName mafplit/ 7species.final.maf

##phyloFit提取4d位点构建,获得nonconserved.mod
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

#根据每条拆分的染色体获得的4d位点文件生成对应染色体的nonconserved.omd
for i in Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC; 
do 
    phyloFit --tree "(oryza_sativa,((((sonneratia_alba,eucalyptus_grandis),arabidopsis_thaliana),(populus_trichocarpa,rhizophora_apiculata)),(sesamum_indicum,avicennia_marina)))" --msa-format SS --out-root $i.nonconserved-4d.mod $i.4d-sites.ss;
done;









#报错
命令行：phastCons --target-coverage 0.25 --expected-length 12 --rho 0.4 --estimate-rho "(oryza_sativa,((((sonneratia_alba,eucalyptus_grandis),arabidopsis_thaliana),(populus_trichocarpa,rhizophora_apiculata)),(sesamum_indicum,avicennia_marina)))" --msa-format MAF Chr1.maf 4dss/Chr1.nonconserved-4d.mod.mod --no-post-probs

Reading tree model from 4dss/Chr1.nonconserved-4d.mod.mod...
Reading alignment from Chr1.maf...
Creating 'conserved' and 'nonconserved' states in HMM...
Computing emission probs (state 0, cat 0, mod 0)...
Computing emission probs (state 1, cat 1, mod 1)...
Finding MLE for (rho)...
ERROR mm_exp_real: got P->size=4, Q->sizse=4, t=-6884.980201


#估计rho来计算myrho.cons.mod和myrho.noncons.mod，随后生成cons.bed
1.计算myrho.cons.mod和myrho.noncons.mod
phastCons --target-coverage 0.25 --expected-length 12 --rho 0.1 --estimate-rho myrho --msa-format MAF Chr1.maf 4dss/all.nonconserved-4d.mod --no-post-probs

2.计算cons elements
phastCons --target-coverage 0.25 --expected-length 12 --msa-format MAF Chr1.maf --most-conserved phastCons.Chr1.bed  myrho.cons.mod,myrho.noncons.mod > phastCons.Chr1.wig

........element.bed结果为空........

#根据 http://compgen.cshl.edu/phast/phastCons-HOWTO.html 中的4.1部分教程重新分析
1.Split the alignments into small fragments

 mkdir -p CHUNKS            # put fragments here
        for file in chr*.maf ; do
            root=$(basename ${file} .maf)
            msa_split ${file} --in-format MAF --refseq ${root}.fa \
                --windows 1000000,0 --out-root CHUNKS/${root} --out-format SS \
                --min-informative 1000 --between-blocks 5000 
        done

2.Estimate parameters for each fragment.

mkdir -p TREES     # put estimated tree models here
rm -f TREES/*      # in case old versions left over
for file in CHUNKS/Chr*.*.ss ; do 
     root=$(basename ${file} .ss) 
     phastCons --target-coverage 0.125 --expected-length 20 \
     --gc 0.4 --estimate-trees TREES/${root} \
     ${file} all.nonconserved-4d.mod --no-post-probs
done

3.Combine the separately estimated parameters by phyloBoot

ls TREES/*.cons.mod > cons.txt
phyloBoot --read-mods '*cons.txt' --output-average ave.cons.mod 
ls TREES/*.noncons.mod > noncons.txt
phyloBoot --read-mods '*noncons.txt' --output-average ave.noncons.mod

4.Predict conserved elements and conservation scores globally using the combined estimates

mkdir -p ELEMENTS SCORES
rm -f ELEMENTS/* SCORES/*
for file in CHUNKS/Chr*.*.ss ; 
do 
    root=$(basename ${file} .ss) 
    phastCons --target-coverage 0.125 --expected-length 20 \
    --most-conserved ELEMENTS/${root}.bed --score \
    ${file} ave.cons.mod,ave.noncons.mod > SCORES/$root.wig
done

......输出结果中ELEMENTS/${root}.bed 全部为空.....