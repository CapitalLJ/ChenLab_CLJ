#!/usr/bin bash

adapter1=AGATCGGAAGAGCACA 
adapter2=AGATCGGAAGAGCGT


#   该脚本的输入文件格式为测序文件,对测序文件进行质控、去接头、提取UMI、bowtie2比对、UMI提取、bam文件的处理（sort、index、mpideup）
# VarScan提取SNP。
#   使用该脚本新建一个文件夹，该文件夹目率下还需要一个VarScan脚本，配置文件脚本、关于配置文件的说明在脚本中会明确说明，用于对文件处理中的一些参数
# 的调整，包括接头序列、UMI的提取、VarScan的一些参数。
# eg：  bash -f gene.fa -1 test-R1.fastq -2 test-R2.fastq   (双端测序)
#       bash -f geng.fa -U test.fastq  (单端测序)

#     建议新建一个文件夹，将需要处理的测序文件移至该目录路径下。会生成一个fastqc文件夹，存放质控数据；过程文件会存放在log文件夹中。UMI_ectract文件夹存放比对前的fastq文件。
# align文件中包括比对后的sam文件、bam文件、bai文件等.UMI_dedup文件夹中的为umi去重后的文件，result文件为所需的vcf文件。

mkdir -p fastqc log UMI_ectract UMI_dedup result adapter blast


while getopts "f:1:2:U:" opt; do
  case $opt in
    f)
      gene=$OPTARG
      ;;
    1)
      file1=$OPTARG
      ;;
    2)
      file2=$OPTARG
      ;;
    U)
      file3=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [[ -n "$file1" && -n "$file2" ]]; then
  pair=1
elif [[ -n "$file3" ]]; then
  pari=0
else
  echo "Usage: $0 [-1 file1 -2 file2] [-U file3]" >&2
  exit 1
fi

name1="${file1%%.*}"   # 去掉文件扩展名
name2="${file2%%.*}"   # 去掉文件扩展名
name3="${file3%%.*}"   # 去掉文件扩展名
genome="${gene%%.*}"

#######################################################################################################
# 质控

mkdir -p  fastqc/primary
if (($pair));then
  fastqc -t 2 -o fastqc/primary $file1 $file2 >/dev/null 2>&1
else
  fastqc -o fastqc/primary $file3 >/dev/null 2>&1
fi







########################################################################################################
#  去接头
# echo "<<<start adapter>>>"

mkdir -p adapter
if (($pair));then
  cutadapt -a $adapter1 -A $adapter2 --minimum-length 40 --overlap 4 \
    -o adapter/$name1 -p adapter/$name2 $file1 $file2 > log/${name1}_adapter.log 2>&1
else
  cutadapt -a $adapter1  --minimum-length 40 --overlap 4 \
    -o adapter/$name3 $file3 > log/${name3}_adapter.log 2>&1
fi

echo "<<<adapter done>>>"

########################################################################################################
# UMI提取
# echo "<<<start UMI_ectract>>>"

mkdir -p UMI_ectract
if (($pair));then
  umi_tools extract --extract-method=regex -I adapter/$name1 --read2-in=adapter/$name2 \
    --bc-pattern="^(?P<umi_1>.{5})(?P<discard_1>.{2}).*" \
    --bc-pattern2="^(?P<umi_2>.{5})(?P<discard_2>.{2}).*" \
    -S ./UMI_ectract/${name1}_umi.gz \
    --read2-out=./UMI_ectract/${name2}_umi.gz > log/${name1}_umi-extract.log 2>&1
else
  umi_tools extract --extract-method=regex --stdin=adapter/$name3 \
    --bc-pattern="^(?P<umi_1>.{5})(?P<discard_1>.{2}).*" \
    --stdout=./UMI_ectract/${name3}_umi.gz > log/${name3}_umi-extract.log 2>&1
fi

echo "<<<UMI_ectract done>>>"

########################################################################################################
# 比对前质控

mkdir -p fastqc/UMI_ectract

if (($pair));then
  fastqc -t 2 -o fastqc/UMI_ectract ./UMI_ectract/${name1}_umi.gz ./UMI_ectract/${name2}_umi.gz >/dev/null 2>&1
else
  fastqc -o fastqc/UMI_ectract ./UMI_ectract/${name3}_umi.gz >/dev/null 2>&1
fi




########################################################################################################
# 构建索引
# echo "<<<start build index>>>"

mkdir -p index
cd index

bowtie2-build ../$gene $genome >/dev/null 2>&1
cd ..

echo "<<<build index done>>>"

########################################################################################################
#比对
# echo "<<<start blast>>>"

mkdir -p blast
if (($pair));then
  bowtie2 -p 24 -x ./index/$genome -1 ./UMI_ectract/${name1}_umi.gz \
   -2 ./UMI_ectract/${name2}_umi.gz -S ./blast/${name1}.sam > log/${name1}_blast.log 2>&1
else
  bowtie2 -p 24 -x ./index/$genome -U ./UMI_ectract/${name3}_umi.gz \
    -S ./blast/${name3}.sam > log/${name3}_blast.log 2>&1
fi

echo "<<<blast done>>>"

########################################################################################################
#UMI去重
# echo "<<<start UMI_dedup>>>"

mkdir -p UMI_dedup
cd ./blast
if (($pair));then
  samtools sort ${name1}.sam -o ${name1}.bam
  samtools index ${name1}.bam
else
  samtools sort ${name3}.sam -o ${name3}.bam
  samtools index ${name3}.bam
fi

if(($pair));then
  umi_tools dedup -I ${name1}.bam --paired --output-stats=../UMI_dedup/${name1} \
    --log=../log/${name1}_UMI_dedup.log -S ../UMI_dedup/${name1}_dedup.bam
else
  umi_tools dedup -I ${name3}.bam  --output-stats=../UMI_dedup/${name3} \
    --log=../log/${name1}_UMI_dedup.log -S ../UMI_dedup/${name3}_dedup.bam  
fi

cd ..

echo "<<<UMI_dedup done>>>"


########################################################################################################
#突变筛选
# echo "<<<start SNP_intarct>>>"

if(($pair));then
  samtools mpileup -q 1 -d 30000 -f $gene ./UMI_dedup/${name1}_dedup.bam \
    1>./UMI_dedup/${name1}_mpileup 2>./log/${name1}_mpileup.log
else
  samtools mpileup -q 1 -d 30000 -f $gene ./UMI_dedup/${name3}_dedup.bam \
    1>./UMI_dedup/${name3}_mpileup 2>./log/${name3}_mpileup.log 
fi


if(($pair));then
  java -jar ./VarScan.v2.4.0.jar pileup2snp ./UMI_dedup/${name1}_mpileup --output-vcf \
    1>./result/${name1}_0.01.snp.vcf 2>log/${name1}_VarScan.log 
else
  java -jar ./VarScan.v2.4.0.jar pileup2snp ./UMI_dedup/${name3}_mpileup --output-vcf \
    1>./result/${name3}_0.01.snp.vcf 2>log/${name3}_VarScan.log 
fi
echo "<<<SNP_intarct done>>>"


























































