#lastz做全基因组比对
nohup lastz /public1/users/zengws/genomedb/arabidopsis_thaliana/Athaliana_447_TAIR10.fa[multiple] /public1/users/zengws/genomedb/rhizophora_apiculata/RA_final.fa[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=maf > at_vs_ra.maf &

nohup lastz /public1/users/zengws/genomedb/arabidopsis_thaliana/Athaliana_447_TAIR10.fa[multiple] /public1/users/zengws/genomedb/sonneratia_alba/SA_final.fa[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=maf > at_vs_sa.maf &

nohup lastz /public1/users/zengws/genomedb/arabidopsis_thaliana/Athaliana_447_TAIR10.fa[multiple] /public1/users/zengws/genomedb/populus_trichocarpa/ptrichocarpa_v4.fa[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=maf > at_vs_pot.maf &

nohup lastz /public1/users/zengws/genomedb/arabidopsis_thaliana/Athaliana_447_TAIR10.fa[multiple] /public1/users/zengws/genomedb/e_grandis/egrandis_v2.fa[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=maf > at_vs_eug.maf &

nohup lastz /public1/users/zengws/genomedb/arabidopsis_thaliana/Athaliana_447_TAIR10.fa[multiple] /public1/users/zengws/genomedb/se_indicum/se.indicum.dna.fa[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=maf > at_vs_sei.maf & 

nohup lastz /public1/users/zengws/genomedb/arabidopsis_thaliana/Athaliana_447_TAIR10.fa[multiple] /public1/users/zengws/genomedb/oryza_sativa/osativa_v7.fa[multiple] --notransition --ambiguous=iupac --step=20 --nogapped --format=maf > at_vs_ors.maf &

#last包下的last-dotplot 绘制可视乎比对点图
last-dotplot at_vs_ors.maf at_vs_ors.png

#