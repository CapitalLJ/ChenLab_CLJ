# nwr 软件说明

nwr is a command line tool for NCBI taxonomy, assembly reports and Newick files, written in Rust.
---nwr是一个用Rust编写的用于NCBI分类法、装配报告和Newick文件的命令行工具。

### Install
Current release: 0.5.10

```shell
brew install wang-q/tap/nwr

```

### SYNOPSIS

```shell
$ nwr help
`nwr` is a command line tool for NCBI taxonomy, assembly reports and Newick files

Usage: nwr [COMMAND]

Commands:
  append    Append fields of higher ranks to a TSV file
  ardb      Init the assembly database
  assembly  Prepare ASSEMBLY materials
  download  Download the latest releases of `taxdump` and assembly reports
  info      Information of Taxonomy ID(s) or scientific name(s)
  lineage   Output the lineage of the term
  member    List members (of certain ranks) under ancestral term(s)
  restrict  Restrict taxonomy terms to ancestral descendants
  txdb      Init the taxonomy database
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

```
#### 准备工作
在使用nwr时，首先使用 `nwr download` 命令来获得NCBI上的相关信息，再使用`nwr txdb`。

```shell
download命令下载的相关文件，在`~/.nwr`中可以查看

# taxdump
wget -N -P ~/.nwr https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget -N -P ~/.nwr https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5

# assembly reports
wget -N -P ~/.nwr https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
wget -N -P ~/.nwr https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

```
txbd命令生成的相关数据库文件，在`~/.nwr`中可以查看，生成的.dmp的说明在[readme.txt](./readme.txt)文件中可以查看。


#### 使用

```shell
nwr info 4932
# 提取ID是4932的物种的信息，注意需要先完成nwr准备工作建立数据库。因为NCBI的信息会更新，所以建议在一个每个project不要频繁使用`nwr download`和`nwr txdb`导致本地数据库版本混乱。下面是上述命令输出结果：

# Saccharomyces cerevisiae - species
# ----------------------------------
# NCBI Taxonomy ID: 4932
# Same as:
# * Candida robusta
# * Mycoderma cerevisiae
# * Saccharomyces capensis
# * Saccharomyces italicus
# * Saccharomyces oviformis
# * Saccharomyces uvarum var. melibiosus
# Commonly named baker's yeast.
# Also known as:
# * brewer's yeast
# * S. cerevisiae
# First description:
# * Mycoderma cerevisiae Desm., 1827
# * Saccharomyces cerevisiae (Desm.) Meyen, 1838
# Part of the Plants and Fungi.

nwr lineage 4932
# 提取ID是4932的物种的系统发育树上各term的信息


# no rank root    1
# no rank cellular organisms      131567
# superkingdom    Eukaryota       2759
# clade   Opisthokonta    33154
# kingdom Fungi   4751
# subkingdom      Dikarya 451864
# phylum  Ascomycota      4890
# clade   saccharomyceta  716545
# subphylum       Saccharomycotina        147537
# class   Saccharomycetes 4891
# order   Saccharomycetales       4892
# family  Saccharomycetaceae      4893
# genus   Saccharomyces   4930
# species Saccharomyces cerevisiae        4932

nwr restrict "Vertebrata" -c 2 -f tests/nwr/taxon.tsv

#筛选tsv文件中包含限定term"Vertebrata"的物种信息
# -c, --column <column> ID所在的列，从1开始[默认：1] 。
# -e, --exclude 排除符合条件的行

nwr member "Homo"

# List members (of certain ranks) under ancestral term(s)
# 输出文件为四列：
# 第一列为ID tax_id
# 第二列是名称 sci_name
# 第三列是rank，等级（中还是亚种）
# 第四列是划分 division

1425170 Homo heidelbergensis  species Primates
63221   Homo sapiens neanderthalensis   subspecies      Primates



nwr append tests/nwr/taxon.tsv -c 2 -r species -r family --id

# 给已有的tsv文件增加描述分类信息

# =>.<=llj:~$ nwr append 11111111.tsv -c 2 -r species -r family --id
# #sci_name       tax_id  species species_id      family  family_id
# Human   9606    Homo sapiens    9606    Hominidae       9604
# Yeast   559292  Saccharomyces cerevisiae        4932    Saccharomycetaceae      4893
# Actinophage JHJ-1       12347   Actinophage JHJ-1       12347   NA      0


```
