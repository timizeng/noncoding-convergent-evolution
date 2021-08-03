#前处理，将3个红树物种蛋白fa序列名前都加上物种标签,#格式sed 's/原字符串/新字符串/' 文件,  sed 's/原字符串/新字符串/g' 文件, 无g代表替换第一个，有g代表替换所有
sed -i 's/>/>Av/g' av.marina.faa   
sed -i 's/>/>So/g' so.alba.faa
sed -i 's/>/>Ra/g' rh.apiculata.faa

#运行orthofinder,寻找同源group
nohup orthofinder -f protein/ -a 12 -t 12 &

#筛选orthfinder 根据输出文件orthogroups.genecount.txt，筛选7个物种均有拷贝（即无 0 copy），且至少有一个物种包含一个copy的group作为备选
#R script
a <- read.csv(file = 'Orthogroups.GeneCount.tsv',header = T,sep = '\t')
zero <- as.data.frame(which(a==0,arr.ind = T))
zero_u <- unique(zero$row)
b <- a[-zero_u,]
one <- as.data.frame(which(b==1,arr.ind = T))
one_u <- unique((one$row))
genecount_one <- b[one_u,]
write.table(genecount_one,file = 'orthogroupone.csv',quote = F,sep = '\t')

#Linux端，制作待计算的filname文件
awk '{print $2}' orthogroupone.csv > countgroup.txt
sed -i '1d' countgroup.txt #删掉第一行的表头
sed -i 's/$/&.fa/g' countgroup.txt #每一个Group后加 .fa 即代表该蛋白文件



#使用awk脚本删除所有序列后的*号
sed -i 's/*//g' input.fa
#使用muscle将每一orthogroup蛋白对齐,随后使用clustalo得到相似度矩阵
muscle -in input.fa -clwstrictout aligned.clw
clustalo -i aligned.clw --percent-id --distmat-out=pim.csv --full --force
sed -i '1d' pim.csv

#使用相似度矩阵进行求和，和即为该序列与其他所有序列的相似度总分,R脚本 orthlist.r
#!/miniconda3/envs/R/bin/Rscript Rscript
library(dplyr)
pim <- read.table(file = 'pim.csv',header = F)
colnames(pim) <- c('name',pim$V1)
species <- unique(substr(pim$name,1,2))
one_group <- c()
newpim <- c()
list <- c()

for (i in 1:length(species)) {
  subone <- pim[grep(pattern = species[i],pim[,1]),]
  if (nrow(subone)==1) {
    newpim <- rbind(subone,newpim)
    one_group <- c(one_group,species[i])
  }
}

rownames(pim) <- pim$name
rownames(newpim) <- newpim$name

#枚举出所有可能组合
species <- setdiff(species, one_group)
species <- gtools::permutations(n=length(species),r=length(species),v=species)

for (k in 1:nrow(species)) {
  multi_species <- species[k,]
  one_pim <- newpim
  for (j in 1:length(multi_species)) {
    subgroup <- pim[grep(pattern = multi_species[j],pim[,1]),]
    target <- as.data.frame(apply(pim[subgroup$name,newpim$name,drop=F],1,sum)) #drop参数保留行列名
    target_top <- target %>% top_n(n=1) %>% rownames()
    top_one <- target_top[1]
    top_one <- pim[grep(pattern = top_one,pim[,1]),]
    one_pim <- top_one %>% full_join(one_pim)
  }
  list <- rbind.data.frame(list,c(one_pim$name,sum(one_pim[,one_pim$name])))
}

final_list <- list[sample(nrow(list),1),]
write.table(final_list,file = 'out.csv',append = T,row.names = T,col.names = F)

#linux端整合脚本 nature.sh
for filename in $(cat countgroup.txt);
do 
   #使用awk脚本将所有序列后的*号
    sed -i 's/*//g' $filename
    #使用muscle将每一orthogroup蛋白对齐,随后使用clustalo得到相似度矩阵
    muscle -in $filename -clwstrictout aligned.clw
    clustalo -i aligned.clw --percent-id --distmat-out=pim.csv --threads=12 --full --force
    sed -i '1d' pim.csv
    Rscript --save orthlist.r
done;

#在orthofinder文件夹下运行脚本,conda的R 4.x环境
conda activate R
bash nature.sh

#输出out.csv进行排序,打开excel→转置矩阵→全选中数据→查看宏→新建宏→复制代码保存执行→转置矩阵
sub colsort()
Set ss = Selection 
For i = 1 To ss.Columns.Count 
ss.Columns(i).Sort Key1:=ss.Columns(i), Order1:=xlAscending, Header:=xlNo, _ 
  Orientation:=xlTopToBottom 
Next 
End Sub











