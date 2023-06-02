



















```shell

# =============================================
# 统计每一个突变的样本个数，下面以X染色体为例
cd ~/llj/human-RNA-structure/SNP
mkdir -p tmp/chrX

parallel -j 24 "
        perl snp_intrac.pl chrX/{1} > chrX/{1}.tsv
" ::: $(ls chrX)


parallel -j 24 "
        cut -f 1-10 {1} > tmp/{1}
" ::: $(ls chrX/*.tsv)


cd tmp
parallel -j 24 "
    perl 1.pl  chrX/gene.list {1} ../output/{1}
" ::: $(ls chrX/*.tsv)



```



```shell
# ============================================================
# 从COSMIC上下载的癌症数据库的处理,提取特定部位癌症的突变数据，并分析同义突变所占的比例。

原文件中的第一、



ath-miR159a uuuggauugaagggagcucua 
hsa-miR-4789-5p GUAUACACCUGAUAUGUGUAUG 
hsa-miR-1200 cuccugagccauucugagccuc 
hsa-miR-146a-3p CCUCUGAAAUUCAGUUCUUCAG 






