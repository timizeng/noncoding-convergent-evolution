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
    echo "nohup lastz-1.04.00 arabidopsis_thaliana.mask.fa[multiple] ${i}[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=axt > ${base}.axt &" >> lastz.sh
done
chmod 755 lastz.sh
./lastz.sh

#axt to psl
for i in *.axt;
do 
    base=$(basename ${i} .axt)
    axtToPsl ${i} /data/users/zengws/conelement/con_ele/chromsize/arabidopsis_thaliana.chromsizes /data/users/zengws/conelement/con_ele/chromsize/${base}.chromsizes ${base}.psl
done

#axtChain
for i in *.psl; 
do 
    base=$(basename ${i} .psl)
    axtChain -linearGap=loose -psl ${i} -faT /data/users/zengws/conelement/con_ele/trfgenome/arabidopsis_thaliana.mask.fa -faQ /data/users/zengws/conelement/con_ele/trfgenome/${base}.mask.fa ${base}.chain
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
	chainPreNet ${i} /data/users/zengws/conelement/con_ele/chromsize/arabidopsis_thaliana.chromsizes /data/users/zengws/conelement/con_ele/chromsize/${base}.chromsizes ${base}.prenet
done

#c:chainNet
for i in *.prenet;
do
    base=$(basename ${i} .prenet)
    chainNet ${i} /data/users/zengws/conelement/con_ele/chromsize/arabidopsis_thaliana.chromsizes /data/users/zengws/conelement/con_ele/chromsize/${base}.chromsizes ${base}.at.net ${base}.net
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
    netToAxt ${i} ${base}.prenet ../twobitgenome/arabidopsis_thaliana.2bit ../twobitgenome/${base}.2bit ${base}.net2axt
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
    axtToMaf -tPrefix=arabidopsis_thaliana. -qPrefix=${base}. ${i} ../chromsize/arabidopsis_thaliana.chromsizes ../chromsize/${base}.chromsizes ${base}.arabidopsis_thaliana.sing.maf 
done

#重命名一下，必须要把arabidopsis_thaliana在前
for i in *.sing.maf; do base=$(basename ${i} .arabidopsis_thaliana.sing.maf); mv ${i} arabidopsis_thaliana.${base}.sing.maf; done;

#step8:combine multiple alignment
roast + E=arabidopsis_thaliana "(((arabidopsis_thaliana theobroma_cacao) (citrus_sinensis (swietenia_macrophylla (xylocarpus_granatum (xylocarpus_rumphii xylocarpus_moluccensis))))) (populus_trichocarpa (carallia_pectinifolia ((bruguiera_gymnorhiza bruguiera_sexangula) ((ceriops_tagal (kandelia_obovata kandelia_candel)) (rhizophora_mangle (rhizophora_apiculata (rhizophora_mucronata rhizophora_stylosa))))))))" *.*.sing.maf 18species.final.maf

