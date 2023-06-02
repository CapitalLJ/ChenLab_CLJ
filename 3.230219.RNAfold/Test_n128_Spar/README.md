以127个酿酒酵母＋S288c（参考）为例，以奇异酵母作为外类群分析SNP与mRNA二级结构的关系

### 全基因组比对
后续补充
### SNP提取
后续补充


### 二级结构预测

预测突变后的mRNA的二级结构以及参考序列的RNA
```shell
mkdir ~/RNA-sumfold
cd ~/RNA-sumfold
mkdir fold ps pdf test result

# 建议在每一步前用test前缀的文件调试一下，运行结果可以储存在test文件中



```

```shell
cd ~/RNA-sumfold
cat Scer_n128_Spar.intact.lst | while read id 
do
mkdir -p fold/"$id"/S288c
mkdir -p ct/"$id"/S288c
done
# mv fold/$id_*.fold fold/"$id"/

# 对每一个突变预测其二级结构
# 后续补充


cat Scer_n128_Spar.intact.lst | while read id
do
cat n128_Spar.gene.fa | grep -A 1 ${id}$ > $id.tmp.fa
RNAfold -T 30 < $id.tmp.fa > fold/"$id"/S288c/$id.fold
rm $id.tmp.fa
done
# 将作为参考的S288c的mRNA的二级结构预测出来后面进行对比
```


### 二级结构比较
首先比较SNP产生后整个二级结构是否改变

```shell
cat Scer_n128_Spar.intact.lst | while read id
do
    cd fold/$id
    ls *.fold | while read name
    do
        tmpname="${name%.fold}"
        IFS='_' read -r site1 site2 site3 site4 <<<"$tmpname"
        A=0
        if [ "$(sed -n '3p' "$name")" = "$(sed -n '3p' "./S288c/$id.fold")" ]; then
            A=1
        fi
        echo -e "${site1}\t${site2}\t${site3}\t${site4}\t${A}">> ~/RNA-sumfold/result/n127_spar_global_compare.tsv
    done
    cd ../../
done
sed -i "1igene\tpos\tREF\tALT\tcompare" ~/RNA-sumfold/result/n127_spar_global_compare.tsv
# 输出文件格式如下
# gene    pos     REF     ALT     compare 
# YAL007C 143     A       G       0       
# YAL007C 157     C       T       1       
# YAL007C 183     T       C       0       
# compare 0表示结构改变，1表示不变



# 将fold文件转化为ct文件来进行后面的比对

cat Scer_n128_Spar.intact.lst | while read id
do
ls fold/"$id"/*.fold | xargs -P 16 -I {} sh -c 'dot2ct {} ct/'$id'/$(basename {} .fold).ct'
ls fold/"$id"/S288c/*.fold | xargs -I {} sh -c 'dot2ct {} ct/'${id}'/S288c/$(basename {} .fold).ct'
done




cat Scer_n128_Spar.intact.lst | while read id
do
    cd ct/$id
    ls *.ct | while read name  
    do
        tmpname="${name%.ct}"
        IFS='_' read -r site1 site2 site3 site4 <<<"$tmpname"
        B=$((site2+1))
        diff <(sed -n "${B}p" $name) <(sed -n "${B}p" ./S288c/$id.ct) | sed -n '2p;4p'  | sed 's/[[:blank:]]\{1,\}/\t/g' | cut -f 6 > $name.tmp
        site5=$(sed -n "2p" $name.tmp | cut -f 1 )
        site6=$(sed -n "1p" $name.tmp | cut -f 1 )
        echo -e "${site1}\t${site2}\t${site3}\t${site4}\t${site5}\t${site6}">> ~/RNA-sumfold/result/n127_spar_compare.tsv
    done
    rm *.tmp
    cd ../../
done
sed -i "1igene\tpos\tREF\tALT\tbefore\tafter" ~/RNA-sumfold/result/n127_spar_compare.tsv

# 将突变后的RNA二级结构和S288c的进行比对，输出文件格式如下
# gene    pos     REF     ALT     before  after
# YAL007C 143     A       G       0       148
# YAL007C 157     C       T       0       0
# YAL007C 183     T       C       194     172
# YAL007C 254     A       G       409     128
# YAL007C 352     A       G       359     359
# YAL007C 412     A       G       250     252
# YAL007C 449     G       T       128     433
# YAL007C 61      G       A       68      0
# YAL007C 648     G       A       0       615

# 第一列表示基因名称、第二列突变在基因上的位置、第五列表示参考碱基的配对情况、第六列表示突变后的配对情况。

```

### freq信息提取

```shell
cd ~/RNA-sumfold/result
mkdir freq_each freq_10

# 根据同义突变和非同义突变、是否改变二级结构、GC获得和失去来对SNP进行筛选。



# 统计S288c的mrna二级结构茎和环的数量,并计算突变后mrna二级结构茎和环的变化
cat Scer_n128_Spar.intact.lst | while read id
do
    cd ct/"$id"/S288c
    ls *.ct | while read name  
    do
        stem=$(sed 's/[[:blank:]]\{1,\}/\t/g' ${name} | cut -f 6 | grep -v ^0 | wc -l )
        ((stem--))
        loop=$(sed 's/[[:blank:]]\{1,\}/\t/g' ${name} | cut -f 6 | grep -w 0 | wc -l )
        echo -e "${id}\t${stem}\t${loop}"> ./stem-loop.tsv
        sed -i "1igene\tstem\tloop" ./stem-loop.tsv  
        # S288c中rna的茎和环的数量
    done
    cd ../../../
    cd ct/$id
    ls *.ct | while read name  
    do
        tmpname="${name%.ct}"
        IFS='_' read -r site1 site2 site3 site4 <<<"$tmpname"
        site5=$(sed 's/[[:blank:]]\{1,\}/\t/g' ${name} | cut -f 6 | grep -v ^0 | wc -l )
        ((site5--))
        site6=$(sed 's/[[:blank:]]\{1,\}/\t/g' ${name} | cut -f 6 | grep -w 0 | wc -l )
        stem=$(sed 's/[[:blank:]]\{1,\}/\t/g' ./S288c/${id}.ct | cut -f 6 | grep -v ^0 | wc -l )
        ((stem--))
        loop=$(sed 's/[[:blank:]]\{1,\}/\t/g' ./S288c/${id}.ct | cut -f 6 | grep -w 0 | wc -l )
        ((site7=site5-stem))
        ((site8=site6-loop))
        echo -e "${site1}\t${site2}\t${site3}\t${site4}\t${site5}\t${site6}\t${site7}\t${site8}">> ~/RNA-sumfold/result/n127_Spar/n127_spar_stem-loop.tsv
    done
    cd ../../

done
sed -i "1igene\tpos\tREF\tALT\tstem\tloop\tstem-change\tloop-change" ~/RNA-sumfold/result/n127_Spar/n127_spar_stem-loop.tsv

# 按照茎是增加还是减少还是不变来进行分类并计算进化速率

cd ~/RNA-sumfold/result/n127_Spar

sed '1d' n127_spar_stem-loop.tsv > tmp.tsv

awk '$7>0 {print}' tmp.tsv | cut -f 10 | perl 1.pl > freq_each/stem-loop/stemup.tsv
awk '$7>0 {print}' tmp.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/stemup.tsv
awk '$7<0 {print}' tmp.tsv | cut -f 10 | perl 1.pl > freq_each/stem-loop/stemdown.tsv
awk '$7<0 {print}' tmp.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/stemdown.tsv
awk '$7==0 {print}' tmp.tsv | cut -f 10 | perl 1.pl > freq_each/stem-loop/stem=.tsv
awk '$7==0 {print}' tmp.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/stem=.tsv

# 从同义突变中按照茎是增加还是减少来进行分类并计算进化速率
cd ~/RNA-sumfold/result/n127_Spar
cut -f 1,2 ./freq_each/sty-nsy/data_syn.tsv > syn.list
grep -w -f syn.list tmp.tsv > n127_spar_syn_stem-loop.tsv

awk '$7>0 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 1.pl > freq_each/stem-loop/syn-stemup.tsv
awk '$7>0 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/syn-stemup.tsv
awk '$7<0 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 1.pl > freq_each/stem-loop/syn-stemdown.tsv
awk '$7<0 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/syn-stemdown.tsv
awk '$7==0 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 1.pl > freq_each/stem-loop/syn-stem=.tsv
awk '$7==0 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/syn-stem=.tsv

# 按照茎增加

awk '$9>0.001 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/0.001-stemup.tsv
awk '$9>0.002 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/0.002-stemup.tsv
awk '$9>0.003 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/0.003-stemup.tsv
awk '$9>0.004 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/0.004-stemup.tsv
awk '$9>0.005 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/0.005-stemup.tsv
awk '$9>0.01 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/0.01-stemup.tsv
awk '$9>0.02 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/0.02-stemup.tsv
awk '$9>0.03 {print}' n127_spar_syn_stem-loop.tsv | cut -f 10 | perl 2.pl > freq_10/stem-loop/0.03-stemup.tsv




for i in YBR145W YBR196C YDR050C YGL253W YGL256W YGR254W YKL060C YLR044C YLR134W YMR083W YMR205C YOR347C YOR374W YPL061W;
do
sum=$(awk '{sum += $2} END {print sum}' ${i}_10)
awk -v sum="$sum" '{print $2/sum}' ${i}_10 | tr '\n' ' ' >> 111.txt
echo >> 111.txt
done



rm tmp.tsv
```
### 进化速率γ拟合

使用matlab工具运用最小二乘法确定进化速率γ。具体操作见get-γ.m。
```bash
for i in YBR145W YBR196C YDR050C YGL253W YGL256W YGR254W YKL060C YLR044C YLR134W YMR083W YMR205C YOR347C YOR374W YPL061W;
do
cat data_SNPs_PARS_syn_codon.csv | grep ${i} | awk -F',' '{print $9"\t"$4"\t"$10"\t"$11"\t"$13"\t"$17}'
done
