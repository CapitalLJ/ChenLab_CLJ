#!/usr/bin bash

# eg: bash vcf-motif-fa.sh  gene.fa file1 file2
# 该脚本的输入文件格式为vcf,该脚本的目的是为了以对照组为参考，统计实验组的高质量的SNP位点并统计其突变情况并为后续寻找motif生成fa文件
# 过程文件会在工作目录生成一个名叫tmp的文件夹里,结果文件在result文件夹中，一个fa文件、两个突变分布文件
# 这个脚本文件需要另外几个perl文件一起，所需的perl文件在Script的perl文件夹里



# filename1：对照组文件的名称
# filename2:实验组文件的名称
# line：文件的行数

mkdir -p tmp
mkdir -p result
gene=$1
filename1=$2
filename2=$3

name1="${filename1%.*}"   # 去掉文件扩展名
name2="${filename2%.*}"   # 去掉文件扩展名

rm -f ./tmp/$name2.list
rm -f ./tmp/$name1-highlevel
rm -f ./tmp/$name2-highlevel



cut -f 1,2,3,4,5,6,7 $filename1 > ./tmp/$name1
cut -f 1,2,3,4,5,6,7 $filename2 > ./tmp/$name2

#########################################################
# 简并碱基替换 （VarScan需要进行）


for i in $name1 $name2;
do
sed -i 's/A\tR/A\tG/g' ./tmp/${i}
sed -i 's/A\tM/A\tC/g' ./tmp/${i}
sed -i 's/A\tW/A\tT/g' ./tmp/${i}
sed -i 's/G\tR/G\tA/g' ./tmp/${i}
sed -i 's/G\tK/G\tT/g' ./tmp/${i}
sed -i 's/G\tS/G\tC/g' ./tmp/${i}
sed -i 's/C\tY/C\tT/g' ./tmp/${i}
sed -i 's/C\tM/C\tA/g' ./tmp/${i}
sed -i 's/C\tS/C\tG/g' ./tmp/${i}
sed -i 's/T\tY/T\tC/g' ./tmp/${i}
sed -i 's/T\tK/T\tG/g' ./tmp/${i}
sed -i 's/T\tW/T\tA/g' ./tmp/${i}
done


##########################################################3
# 对实验组和突变组的SNP进行质量筛选
line3=`cat ./tmp/$name1 | wc -l `
line4=`cat ./tmp/$name2 | wc -l `

for n in $(seq 2 $line3);do
    Chrom1=$(sed -n "${n}p;${n}q" ./tmp/$name1 | cut -f 1)
    Position1=$(sed -n "${n}p;${n}q" ./tmp/$name1 | cut -f 2)
    Ref1=$(sed -n "${n}p;${n}q" ./tmp/$name1 | cut -f 3)
    Cons1=$(sed -n "${n}p;${n}q" ./tmp/$name1 | cut -f 4)
    Reads1=$(sed -n "${n}p;${n}q" ./tmp/$name1 | cut -f 5)
    Reads2=$(sed -n "${n}p;${n}q" ./tmp/$name1 | cut -f 6)
    Reads=$(echo "$Reads1+$Reads2" | bc)
    VarFreq=$(sed -n "${n}p;${n}q" ./tmp/$name1 | cut -f 7)
    VarFreq=$(echo $VarFreq | sed 's/%//g' | awk '{print $1/100}')
    if [ $Reads -ge 10 ];then
        if  (($(echo "$VarFreq > 0.1" | bc))) ;then
            echo -e "$Chrom1\t$Position1\t$Ref1\t$Cons1\t$Reads\t$VarFreq" >> ./tmp/$name1-highlevel
        fi
    elif [ $Reads -ge 8 ];then
        if  (($(echo "$VarFreq > 0.3" | bc))) ;then
            echo -e "$Chrom1\t$Position1\t$Ref1\t$Cons1\t$Reads\t$VarFreq" >> ./tmp/$name1-highlevel
        fi
    fi
done

for n in $(seq 2 $line4);do
    Chrom1=$(sed -n "${n}p;${n}q" ./tmp/$name2 | cut -f 1)
    Position1=$(sed -n "${n}p;${n}q" ./tmp/$name2 | cut -f 2)
    Ref1=$(sed -n "${n}p;${n}q" ./tmp/$name2 | cut -f 3)
    Cons1=$(sed -n "${n}p;${n}q" ./tmp/$name2 | cut -f 4)
    Reads1=$(sed -n "${n}p;${n}q" ./tmp/$name2 | cut -f 5)
    Reads2=$(sed -n "${n}p;${n}q" ./tmp/$name2 | cut -f 6)
    Reads=$(echo "$Reads1+$Reads2" | bc)
    VarFreq=$(sed -n "${n}p;${n}q" ./tmp/$name2 | cut -f 7)
    VarFreq=$(echo $VarFreq | sed 's/%//g' | awk '{print $1/100}')
    if [ $Reads -ge 10 ];then
        if  (($(echo "$VarFreq > 0.1" | bc))) ;then
            echo -e "$Chrom1\t$Position1\t$Ref1\t$Cons1\t$Reads\t$VarFreq" >> ./tmp/$name2-highlevel
        fi
    elif [ $Reads -ge 8 ];then
        if  (($(echo "$VarFreq > 0.3" | bc))) ;then
            echo -e "$Chrom1\t$Position1\t$Ref1\t$Cons1\t$Reads\t$VarFreq" >> ./tmp/$name2-highlevel
        fi
    fi
done

# 提取对照组和实验组的相同突变位点和不同突变位点，以后后续分析
awk 'NR==FNR {a[$1$2$3$4]=$0; next} ($1$2$3$4 in a) {print a[$1$2$3$4]}' ./tmp/$name1-highlevel ./tmp/$name2-highlevel > ./tmp/$name1-commom
awk 'NR==FNR {a[$1$2$3$4]=$0; next} ($1$2$3$4 in a) {print }' ./tmp/$name1-highlevel ./tmp/$name2-highlevel > ./tmp/$name2-commom

cat ./tmp/$name1-highlevel | grep -v -f ./tmp/$name1-commom > ./tmp/$name1-diff
cat ./tmp/$name2-highlevel | grep -v -f ./tmp/$name2-commom > ./tmp/$name2-diff

line1=`cat ./tmp/$name1-commom | wc -l `

line2=`cat ./tmp/$name2-diff | wc -l `

# 相同的位点实验组和对照组做差处理
for n in $(seq 1 $line1);do
Chrom1=$(sed -n "${n}p;${n}q" ./tmp/$name1-commom | cut -f 1)
Position1=$(sed -n "${n}p;${n}q" ./tmp/$name1-commom | cut -f 2)
Ref1=$(sed -n "${n}p;${n}q" ./tmp/$name1-commom | cut -f 3)
Cons1=$(sed -n "${n}p;${n}q" ./tmp/$name1-commom | cut -f 4)
Reads=$(sed -n "${n}p;${n}q" ./tmp/$name2-commom | cut -f 5)


VarFreq1=$(sed -n "${n}p;${n}q" ./tmp/$name1-commom | cut -f 6)
VarFreq1=$(echo $VarFreq1 | sed 's/%//g' | awk '{print $1/100}')

VarFreq2=$(sed -n "${n}p;${n}q" ./tmp/$name2-commom | cut -f 6)
VarFreq2=$(echo $VarFreq2 | sed 's/%//g' | awk '{print $1/100}')

VarFreq=$(echo "scale=2; $VarFreq2-$VarFreq1" | bc)
if [ $Reads -ge 10 ];then
    if  (($(echo "$VarFreq > 0.1" | bc))) ;then
        echo -e "$Chrom1\t$Position1\t$Ref1\t$Cons1\t$Reads\t$VarFreq" >> ./tmp/$name2.list
    fi
elif [ $Reads -ge 8 ];then
    if  (($(echo "$VarFreq > 0.3" | bc))) ;then
        echo -e "$Chrom1\t$Position1\t$Ref1\t$Cons1\t$Reads\t$VarFreq" >> ./tmp/$name2.list
    fi
fi
done


# 实验组特有的突变进行筛选
# for n in $(seq 1 $line2);do
#     Chrom1=$(sed -n "${n}p;${n}q" ./tmp/$name2-diff | cut -f 1)
#     Position1=$(sed -n "${n}p;${n}q" ./tmp/$name2-diff | cut -f 2)
#     Ref1=$(sed -n "${n}p;${n}q" ./tmp/$name2-diff | cut -f 3)
#     Cons1=$(sed -n "${n}p;${n}q" ./tmp/$name2-diff | cut -f 4)
#     Reads=$(sed -n "${n}p;${n}q" ./tmp/$name2-diff | cut -f 5)
#     VarFreq=$(sed -n "${n}p;${n}q" ./tmp/$name2-diff | cut -f 6)
#     VarFreq=$(echo $VarFreq | sed 's/%//g' | awk '{print $1/100}')
#     if [ $Reads -ge 10 ];then
#         if  (($(echo "$VarFreq > 0.1" | bc))) ;then
#             echo -e "$Chrom1\t$Position1\t$Ref1\t$Cons1\t$Reads\t$VarFreq" >> ./tmp/$name2.list
#         fi
#     elif [ $Reads -ge 8 ];then
#         if  (($(echo "$VarFreq > 0.3" | bc))) ;then
#             echo -e "$Chrom1\t$Position1\t$Ref1\t$Cons1\t$Reads\t$VarFreq" >> ./tmp/$name2.list
#         fi
#     fi
# done
cat ./tmp/$name2-diff >> ./tmp/$name2.list



awk '{print $1"\t"($2-50)"\t"($2+50)"\t+"}' ./tmp/$name2.list > ./tmp/$name2.fa.list

perl pick_seq_from_fasta.pl $gene ./tmp/$name2.fa.list > ./result/$name1-$name2.fa

perl 1.pl ./tmp/$name1-highlevel > ./result/$name1-ferq10
perl 1.pl ./tmp/$name2-highlevel > ./result/$name2-ferq10

