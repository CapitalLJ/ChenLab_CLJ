# SNP_pipeline

### 前期准备

创建文件夹以存放不同种类的数据

```shell
cd ~/ecoli_snp
mkdir genome sequence output Script
# sequence 测序数据 genome 菌株参考基因组  
```

### 过滤质控

1、质量评估

```shell
cd ~/ecoli_snp/sequence
mkdir -p ../output/fastqc
fastqc -t 6 -o ../output/fastqc *.gz
cd ../output/fastqc
multiqc

```

2、去接头

NJU8108 NJU8109 NJU8110 NJU8111四个样本为基因组测序
NJU8112 NJU8113 NJU8114 NJU8115四个样本为转录组测序




```shell
# 基因组分析
cd ~/ecoli_snp/sequence
mkdir -p ../output/adapter

for i in NJU8108 NJU8109 NJU8110 NJU8111; do
  cutadapt -a AGATCGGAAGAGCACA -A AGATCGGAAGAGCGT \
    --minimum-length 40 --overlap 4 \
    -o ../output/adapter/${i}_R1.fq.gz -p ../output/adapter/${i}_R2.fq.gz ${i}/${i}_R1.fq.gz ${i}/${i}_R2.fq.gz
done >> ../output/adapter/adapter2320227.log


mkdir -p ../output/fastqc/adapter
cd ~/ecoli_snp/output/adapter 
fastqc -t 6 -o ../output/fastqc/adapter/ *.gz

# 转录组分析
cd ~/ecoli_snp/sequence
mkdir -p ../output/tmr-adapter

for i in NJU8112 NJU8113 NJU8114 NJU8115; do
trimmomatic PE -threads 16 -phred33 \
${i}/${i}_R1.fq.gz ${i}/${i}_R2.fq.gz \
./output/tmr-adapter/${i}_R1_paired.fq.gz ./output/tmr-adapter/${i}_R1_unpaired.fq.gz \
./output/tmr-adapter/${i}_R2_paired.fq.gz ./output/tmr-adapter/${i}_R2_unpaired.fq.gz \
ILLUMINACLIP:/home/linuxbrew/.linuxbrew/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 MINLEN:40 2>&1 >./output/tmr-adapter/${i}.log
done



cd ../output/tmr-adapter
fastqc -t 6 -o ../output/fastqc/tmr-adapter/ *_paired.fq.gz

# 质控后发现5‘端碱基的T含量很高，决定分为剔除和不剔除分别比对一下
for i in  NJU8113 NJU8114 NJU8115; do
  umi_tools extract --extract-method=regex \
  -I ./${i}_R1_paired.fq.gz --read2-in=./${i}_R2_paired.fq.gz \
  --bc-pattern="^(?P<umi_1>.{10}).*$" \
  --bc-pattern2="^(?P<umi_2>.{10}).*$" \
  -S ../umi_ectract_adapter/${i}_R1.fq.gz --read2-out=../umi_ectract_adapter/${i}_R2.fq.gz >>../umi_ectract_adapter/${i}.log
done

umi_tools extract --extract-method=regex \
  -I 1 --read2-in=2 \
  --bc-pattern="^(?P<umi_1>.{10}).*$" \
  --bc-pattern2="^(?P<umi_2>.{10}).*$" \
  -S 11 --read2-out=22 >>1.log


```


3、UMI 提取

```shell
cd ~/ecoli_snp/output/adapter
mkdir ../UMI_ectract

for i in NJU81{09..11}; do
  umi_tools extract --extract-method=regex \
  -I ./${i}_R1.fq.gz --read2-in=./${i}_R2.fq.gz \
  --bc-pattern="^(?P<umi_1>.{5})(?P<discard_1>.{2}).*" \
  --bc-pattern2="^(?P<umi_2>.{5})(?P<discard_2>.{2}).*" \
  -S ../UMI_ectract/${i}_R1.fq.gz --read2-out=../UMI_ectract/${i}_R2.fq.gz >../UMI_ectract/{}.log
done






```

### 构建索引

```shell
# 以大肠杆菌BL21基因组为参考序列 
# https://www.ncbi.nlm.nih.gov/nuccore/NC_012971.2

cd ~/ecoli_snp/genome
mkdir -p index

cd index
bowtie2-build ../BL21.fasta BL21
hisat2-build ../BL21.fasta BL21
# bsub -q serial -n 4 -J "index_build" "bowtie2-build ../BL21.fasta BL21"

mkpir -p A3A-index
mkpir -p APO2-index

# 加入了A3A和APO2的质粒的序列进行比对
```

### 比对

```shell
cd ~/ecoli_snp/output/UMI_ectract
mkdir ../align

bsub -q mpi -n 24 -J "BOWTIE2_alignllj" "bash 1.sh"

bsub -q mpi -n 24 -J "hsiat2_alignllj" "bash hisat2.sh"

# 1.sh
# #!/usr/bin bash
for i in NJU8108 NJU8109 NJU8110 NJU8111;
do
bowtie2 -p 24 -x ../../genome/index/BL21 -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz -S ../align/${i}.sam
done

hisat2.sh
#!/usr/bin bash
for i in NJU8112 NJU8113 NJU8114 NJU8115;
do
hisat2 --threads 20 ../../genome/h-index/BL21 -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz -S ../align/${i}.sam 2>../align/${i}.log
done 


#3371441 reads; of these:
#  3371441 (100.00%) were paired; of these:
#    586849 (17.41%) aligned concordantly 0 times
#    2674951 (79.34%) aligned concordantly exactly 1 time
#    109641 (3.25%) aligned concordantly >1 times
#    ----
#    586849 pairs aligned concordantly 0 times; of these:
#      159297 (27.14%) aligned discordantly 1 time
#    ----
#    427552 pairs aligned 0 times concordantly or discordantly; of these:
#      855104 mates make up the pairs; of these:
#        652438 (76.30%) aligned 0 times
#        180254 (21.08%) aligned exactly 1 time
#        22412 (2.62%) aligned >1 times
#90.32% overall alignment rate

cd ../align
ls
# NJU8108.sam  NJU8109.sam  NJU8110.sam  NJU8111.sam
parallel -j 4 '
  samtools sort {}.sam -o {}.bam
' ::: NJU8112 NJU8113 NJU8114 NJU8115

```

### UMI去重

```shell
cd ~/ecoli_snp/output/align

parallel -j 4 '
  samtools index {}.bam
' ::: NJU8108 NJU8109 NJU8110 NJU8111
# 在UMI去重步骤中索引文件是必须的

for i in NJU81{08..11}; do
  umi_tools dedup -I ${i}.bam --paired --output-stats=../UMI_dedup/${i} \
    --log=../UMI_dedup/${i}.log \
    -S ../UMI_dedup/${i}_dedup.bam
done
# --output-stats 提供一个范围 有关运行的统计信息
# --paired  双端测序选项
# --methed 选择不同的去重方法。 directional（默认选项）：识别连接的 UMI 的聚类（基于汉明距离） 阈值）和 umi A 计数 >= （2* umi B 计数）- 1。

#2023-02-16 16:37:48,428 INFO Reads: Input Reads: 3761213, Read pairs: 3761213, Read 2 unmapped: 144622, Read 1 unmapped: 75672
#2023-02-16 16:37:48,428 INFO Number of reads out: 2520723
#2023-02-16 16:37:48,428 INFO Total number of positions deduplicated: 2498474
#2023-02-16 16:37:48,428 INFO Mean number of unique UMIs per position: 1.01
#2023-02-16 16:37:48,428 INFO Max. number of unique UMIs per position: 6
```

### 突变筛选

```shell
for i in NJU81{12..15};
do
samtools rmdup -S ${i}.bam ${i}-clean.bam
done

```


parallel -j 4 '
  samtools depth {}.bam
' ::: NJU8112 NJU8113 NJU8114 NJU8115






```shell
# 使用VarScan来 call SNP
cd ~/ecoli_snp/output
mkdir VarScan
cd VarScan
mkdir SNP mpileup

cd ~/ecoli_snp/output/UMI_dedup
for i in NJU81{12..15}; do
  samtools mpileup -q 1 -d 30000 -f ../../genome/BL21.fasta \
    ${i}.bam 1>../VarScan/mpileup/${i}.mpileup 2>${i}.log
done





# mpileup文件是VarScan call snp所必需的

cd ../VarScan/mpileup
# 2.sh
# #!/usr/bin bash
# for i in NJU81{08..11}; do
# java -jar ../../../VarScan.v2.4.0.jar pileup2snp ${i}.mpileup --output-vcf 1 >../SNP/${i}_snp.vcf
# done

# bsub -q serial -n 1 -J "call-snp" "bash 2.sh" 

# cd ~/ecoli_snp/output/VarScan/SNP
# ls
# NJU8108_snp.vcf NJU8109_snp.vcf NJU8110_snp.vcf NJU8111_snp.vcf

# --min-var-freq =0.01
for i in NJU81{08..11}; do
java -jar ../../../VarScan.v2.4.0.jar pileup2snp ${i}.mpileup --output-vcf 1 >../SNP/${i}_0.01.snp.vcf
done


for i in NJU81{08..11}; do
java -jar ../../../VarScan.v2.4.0.jar pileup2snp ${i}.mpileup --output-vcf 1 >../SNP/${i}_0.01.snp.vcf
done

parallel -j 4 "
java -jar ../VarScan.v2.4.0.jar pileup2snp {}.mpileup --min-var-freq 0.05 --output-vcf 1 >./{}_0.01.snp.vcf " :::  NJU8112 NJU8113 NJU8114 NJU8115 

parallel -j 4 "
java -jar ../VarScan.v2.4.0.jar pileup2snp {}.mpileup --min-var-freq 0.1 --output-vcf 1 >./{}_0.01.snp.vcf " :::  NJU8112 NJU8113 NJU8114 NJU8115 


# --min-var-freq =0.05
for i in NJU81{08..11}; do
java -jar ../../../VarScan.v2.4.0.jar pileup2snp ${i}.mpileup --min-var-freq 0.05 --output-vcf 1 >../SNP/${i}_0.05.snp.vcf
done

# --min-var-freq =0.1
for i in NJU81{08..11}; do
java -jar ../../../VarScan.v2.4.0.jar pileup2snp ${i}.mpileup --min-var-freq 0.1 --output-vcf 1 >../SNP/${i}_0.1.snp.vcf
done

# 得到SNP的信息以供后续分析

#samtools mpileup -q 1 -d 30000 \
#  --ignore-RG -f DE3.fa \
#  NJU81{08..11}.dedup.bam \
#  >total.mpileup
#varscan mpileup2snp \
#  total.mpileup --output-vcf 1 \
#  --p-value 0.05 --min-var-freq 0.01 \
#  >total.vcf
```

### SNP位点分析

```shell
# 将实验组和对照组相同位点突变信息提取到.comm文件中，实验组新产生的突变信息提取到.diff文件中。
cd ~/ecoli_snp/output/VarScan/SNP
for i in 0.01 0.05 0.1;do
    for x in NJU81{12..15};do
        cut -f 2,3,4 ${x}_${i}-clean.snp.vcf > ${x}_${i}.tmp-clean
        cut -f 2,3,4,7 ${x}_${i}-clean.snp.vcf > ${x}_${i}-clean
    done
    perl 2.pl NJU8112_${i}.tmp-clean NJU8113_${i}.tmp>12_13_${i}both-clean
    perl 2.pl NJU8114_${i}.tmp-clean NJU8115_${i}.tmp>14_15_${i}both-clean

    cut -f 1 12_13_${i}both-clean >12_13_${i}both1-clean
    cut -f 3 12_13_${i}both-clean >12_13_${i}both2-clean
    cut -f 1 14_15_${i}both-clean >14_15_${i}both1-clean
    cut -f 3 14_15_${i}both-clean >14_15_${i}both2-clean

    perl 3.pl 12_13_${i}both1-clean 12_13_${i}both2-clean NJU8112_${i}-clean >../contrast_snp/NJU8112_${i}.comm-clean
    perl 3.pl 12_13_${i}both1-clean 12_13_${i}both2-clean NJU8113_${i}-clean >../contrast_snp/NJU8113_${i}.comm-clean
    perl 3.pl 14_15_${i}both1-clean 14_15_${i}both2-clean NJU8114_${i}-clean >../contrast_snp/NJU8114_${i}.comm-clean
    perl 3.pl 14_15_${i}both1-clean 14_15_${i}both2-clean NJU8115_${i}-clean >../contrast_snp/NJU8115_${i}.comm-clean

    perl 4.pl 12_13_${i}both1-clean 12_13_${i}both2-clean NJU8112_${i}-clean >../contrast_snp/NJU8112_${i}.diff-clean
    perl 4.pl 12_13_${i}both1-clean 12_13_${i}both2-clean NJU8113_${i}-clean >../contrast_snp/NJU8113_${i}.diff-clean
    perl 4.pl 14_15_${i}both1-clean 14_15_${i}both2-clean NJU8114_${i}-clean >../contrast_snp/NJU8114_${i}.diff-clean
    perl 4.pl 14_15_${i}both1-clean 14_15_${i}both2-clean NJU8115_${i}-clean >../contrast_snp/NJU8115_${i}.diff-clean
    rm *.both-clean *.tmp-clean *.both1-clean *.both2-clean
done


for x in NJU81{12..15};do
    cut -f 2,3,4,7 ${x}_0.01-clean-20x.snp.vcf > ${x}_0.01-clean-20x
done


# 将复合碱基替换成突变碱基
for i in NJU81{12..15}
do
sed -i 's/A\tR/A\tG/g' ${i}_0.1-clean
sed -i 's/A\tM/A\tC/g' ${i}_0.1-clean
sed -i 's/A\tW/A\tT/g' ${i}_0.1-clean
sed -i 's/G\tR/G\tA/g' ${i}_0.1-clean
sed -i 's/G\tK/G\tT/g' ${i}_0.1-clean
sed -i 's/G\tS/G\tC/g' ${i}_0.1-clean
sed -i 's/C\tY/C\tT/g' ${i}_0.1-clean
sed -i 's/C\tM/C\tA/g' ${i}_0.1-clean
sed -i 's/C\tS/C\tG/g' ${i}_0.1-clean
sed -i 's/T\tY/T\tC/g' ${i}_0.1-clean
sed -i 's/T\tK/T\tG/g' ${i}_0.1-clean
sed -i 's/T\tW/T\tA/g' ${i}_0.1-clean
done
# 将实验组的突变率减去对照的突变率得到修正后的实验组突变率，将结果整和到.sub文件中 
```









```shell
cd ~/ecoli_snp/output/contrast_snp

# 将突变率按照10％一个区间统计不同区间的突变的信息
for x in 12-13 14-15; do
  for i in {1..10};do
    a=$(echo "scale=1; ($i / 10) - 0.1" | bc)
    b=$(echo "scale=1; ($i / 10)" | bc)
    tsv-filter ${x}-clean.sub -H --ge VarFreq:$a --lt VarFreq:$b > ${x}-clean_$a-$b.tsv
  done
done


# 对于A3A（08-09）选择0.1 0.2 0.3三个不同的最低突变率；AP2选择0.2为最低突变率。
tsv-filter 08_09.sub -H --ge VarFreq:0.1  > 08_09_0.1-1.tsv
tsv-filter 08_09.sub -H --ge VarFreq:0.2  > 08_09_0.2-1.tsv
tsv-filter 08_09.sub -H --ge VarFreq:0.3  > 08_09_0.3-1.tsv
tsv-filter 10_11.sub -H --ge VarFreq:0.2  > 10_11_0.2-1.tsv

tsv-filter 12-13.sub -H --ge VarFreq:0.3  > ../12_13_0.3-1.tsv
tsv-filter 14-15.sub -H --ge VarFreq:0.3  > ../14_15_0.3-1.tsv

tsv-filter 12-13-clean.sub -H --ge VarFreq:0.1  > ../12_13_0.1-1-clean.tsv
tsv-filter 14-15-clean.sub -H --ge VarFreq:0.1  > ../14_15_0.1-1-clean.tsv

tsv-filter 12-13-clean-20x.sub -H --ge VarFreq:0.1  > ../12_13_0.1-1-clean-20x.tsv
tsv-filter 14-15-clean-20x.sub -H --ge VarFreq:0.1  > ../14_15_0.1-1-clean-20x.tsv

# 得到突变位点及其信息后分析其上下游序列是否具有特定的motif，首先提取突变的位点前后50个碱基作为分析序列。

awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 08_09_0.1-1.tsv > 08_09_0.1-1.list
awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 08_09_0.2-1.tsv > 08_09_0.2-1.list
awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 08_09_0.3-1.tsv > 08_09_0.3-1.list
awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 10_11_0.2-1.tsv > 10_11_0.2-1.list

awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 12_13_0.3-1.tsv > 12_13_0.3-1.list
awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 14_15_0.3-1.tsv > 14_15_0.3-1.list

awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 12_13_0.1-1-clean.tsv > 12_13_0.1-1-clean.list
awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 14_15_0.1-1-clean.tsv > 14_15_0.1-1-clean.list

awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 12_13_0.1-1-clean-20x.tsv > 12_13_0.1-1-clean-20x.list
awk '{print "NC_012971.2\t"($1-50)"\t"($1+50)"\t+"}' 14_15_0.1-1-clean-20x.tsv > 14_15_0.1-1-clean-20x.list


perl pick_seq_from_fasta.pl ~/ecoli_snp/genome/BL21.fa 08_09_0.1-1.list > 08_09_0.1-1.fa
perl pick_seq_from_fasta.pl ~/ecoli_snp/genome/BL21.fa 08_09_0.2-1.list > 08_09_0.2-1.fa
perl pick_seq_from_fasta.pl ~/ecoli_snp/genome/BL21.fa 08_09_0.3-1.list > 08_09_0.3-1.fa
perl pick_seq_from_fasta.pl ~/ecoli_snp/genome/BL21.fa 10_11_0.2-1.list > 10_11_0.2-1.fa
```

```Rscript
install.packages('ggseqlogo')
install.packages('seqinr')

library(ggseqlogo)
library(seqinr)
a1 <- read.fasta('D:/08_09_0.1-1.fa')
a2 <- read.fasta('D:/08_09_0.2-1.fa')
a3 <- read.fasta('D:/08_09_0.3-1.fa')
b1 <- read.fasta('D:/10_11_0.2-1.fa')

result <- list()
title <- c("APOBEC3A_0.1", "APOBEC3A_0.2", "APOBEC3A_0.3", "APOBEC2_0.2")
varnames <- c("a1", "a2", "a3", "b1")
for (i in 1:4) {
  temp_str = vector(mode = 'character')
  for (j in 1:length(get(varnames[i]))) {
    temp_str[j] = toupper(c2s(get(varnames[i])[[j]]))
  }
  result[[i]] <- ggseqlogo(temp_str) +
    ggtitle(title[i]) #+
}
gridExtra::grid.arrange(
  result[[1]],
  result[[2]],
  result[[3]],
  result[[4]],
  ncol = 1)

```

2.pl

提取两个文件中的相同行

```perl5
#!/usr/bin/perl
use strict;
use warnings;

open(my $in1, "<", $ARGV[0]);
open(my $in2, "<", $ARGV[1]);

my %compare;
while (<$in1>) {
    chomp;
    $compare{$_} = 1;
}

while (<$in2>) {
    chomp;
    if (exists($compare{$_})) {
        print("$_\n");
    }
}
```

3.pl

按照列表文件提取文件中的目标突变位点

```perl5
#!/usr/bin/perl

my ($fileA, $fileB, $fileC) = @ARGV;

open A, '<', $fileA or die "Unable to open file:$fileA:$!";
while (<A>) {
    chomp;
    push @T, $_;
}
close A;

open B, '<', $fileB or die "Unable to open file:$fileB:$!";
while (<B>) {
    chomp;
    push @B, $_;
}
close B;

$l = $#T;
$b = 0;
# print "$l\n";
open C, '<', $fileC or die "Unable to open file:$fileC:$!";
while (<C>) {
    chomp;
    @N = ($b .. $l);
    for $a (@N) {
        #		print "$a\n";
        $x = @T[$a];
        chomp($x);
        $y = @B[$a];
        chomp($y);
        #		print "$x\n";
        if (/^$x\t[A-Z]\t$y/) {
            print "$_\n";
            $b = $b + 1;
        }
    }
}
close C;
```

4.pl
按照列表文件剔除文件中的目标突变位点

```perl5
#!/usr/bin/perl

my ($fileA, $fileB, $fileC) = @ARGV;

open A, '<', $fileA or die "Unable to open file:$fileA:$!";
while (<A>) {
    chomp;
    push @T, $_;
}
close A;

open B, '<', $fileB or die "Unable to open file:$fileB:$!";
while (<B>) {
    chomp;
    push @B, $_;
}
close B;

$l = $#T;
$b = 0;
# print "$l\n";
open C, '<', $fileC or die "Unable to open file:$fileC:$!";
while (<C>) {
    push @C, $_;
    @N = ($b .. $l);
    for $a (@N) {
        #		print "$a\n";
        $x = @T[$a];
        chomp($x);
        $y = @B[$a];
        chomp($y);
        #		print "$x\n";
        if (/^$x\t[A-Z]\t$y/) {
            pop @C;
            $b = $b + 1;
        }
    }
}
close C;
```

暂时放在这里，，，，，，，，，，，，，

```shell
for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    >&2 echo "==> ${NAME}"
    spanr statop \
        ../blast/S288C.sizes \
        llj-genes.non-overlapped.yml ${NAME}.yml \
        --op intersect --all -o stdout |
        grep -v "^key" |
        perl -nla -F, -e '
            $F[2] == $F[4] and print $F[0];
        ' \
        > llj-${NAME}.intact.lst
done




mkdir -p llj
cat llj_gene_list.csv |
    parallel --colsep ',' --no-run-if-empty --linebuffer -k -j 12 '
        echo {1}
        echo {2} | spanr cover stdin -o llj/{1}.yml
    '

spanr merge llj/*.yml -o llj.merge.yml


for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    >&2 echo "==> ${NAME}"
    mkdir -p llj_${NAME}

    cat llj-${NAME}.intact.lst |
        parallel --no-run-if-empty --linebuffer -k -j 12 "
           fasops slice ${NAME}.fas.gz llj/{}.yml -n S288C -o llj_${NAME}/{}.fas
        "
done

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    >&2 echo "==> ${NAME}"
    mkdir -p llj-SNP_${NAME}

    cat llj-${NAME}.intact.lst |
        parallel --no-run-if-empty --linebuffer -k -j 24 "
            fasops vars --outgroup --nocomplex llj_${NAME}/{}.fas -o stdout |
                sed 's/\$/\t{}/' \
                > llj-SNP_${NAME}/{}.tsv
        "

    #loccation,REF,ALT,mutant_to,freq,occured,gene
    cat llj-SNP_${NAME}/*.tsv |
        tsv-select -f 5,6,7,9,10,8,14 \
        > llj-${NAME}.SNPs.tsv
done

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    tsv-join -z \
        ../vcf/1011Matrix.ext.txt \
        -f ../gene-filter/llj-${NAME}.SNPs.tsv \
        --key-fields 1 \
        --append-fields 2-7 \
        > llj-${NAME}.SNPs.tsv

    cat llj-${NAME}.SNPs.tsv | datamash check

    cat llj-${NAME}.SNPs.tsv |
        perl -nla -F"\t" -e '
            my $loc = $F[0];
            $loc =~ /^(.*):(.*)/;
            my $chr = $1;
            my $pos = $2;
            print qq{$chr\t$pos\t$pos\t$F[1]\t$F[2]};
        ' \
        > llj-${NAME}.upload.tsv
done

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    cat llj-${NAME}.vep.txt |
        perl -nla -F"\t" -e '
            next if /^#/;
            my $loca = $F[1];
            $loca =~ /^(.*)-[0-9]+/;
            my $ID = $1;
            #location,allele,gene,consequence,CDS_position,amino_acids,codons,existing_variation
            print qq{$ID\t$F[2]\t$F[3]\t$F[6]\t$F[8]\t$F[10]\t$F[11]\t$F[12]};
        ' \
    > llj-${NAME}.vep.tsv
done



```
