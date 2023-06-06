###  安装nwr并启动本地数据库

```shell
brew install wang-q/tap/nwr # 0.5.4 or above
brew install sqlite 
# SQLite是一种开源的嵌入式关系型数据库管理系统（DBMS），广泛应用于各种应用程序和平台中。

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank

```
关于nwr的说明，详见[此页面](https://github.com/wang-q/nwr)


# Bacteria
All genomes of Bacteria and Archaea, species by species.Download all genomes and analyze representative strains.

### List all ranks

域（Domain）、界（Kingdom）、门（Phylum）、纲（Class）、目（Order）、科（Family）、属（Genus）、种（Species）

```shell
nwr member Bacteria |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat


# 拉丁名中"sp."是某属中一未知种的意思，"spp."是某属中多个未知种的意思。 "sp."和"spp."是种species的缩写，"sp."是种的单数，"spp."是种的复数。
# tsv-summarize对标签分隔值文件中的字段进行聚合操作。

nwr member Archaea |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat

```

### Species with assemblies
<!-- #region -->
<!-- 提取细菌和古菌的所有属的信息 -->
```shell
mkdir /mnt/Leilingjie/Bacteria/summary
cd /mnt/Leilingjie/Bacteria/summary

nwr member Bacteria Archaea -r genus |
    grep -v -i "Candidatus " |   
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#4470 genus.list
head genus.list.tsv

# genus_id    
# 6       Azorhizobium    genus   Bacteria
# 10      Cellvibrio      genus   Bacteria
# 13      Dictyoglomus    genus   Bacteria
# 16      Methylophilus   genus   Bacteria
```
<!-- #endregion -->



<!-- #region -->
<!-- 下面是四个水平的物种信息提取 -->
```shell


# =============================================================================

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
echo "
    SELECT
        species_id,
        species,
        COUNT(DISTINCT tax_id) AS count -- with strain ID
    FROM ar
    WHERE 1=1
        AND genus_id = ${RANK_ID}
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    GROUP BY species_id
    HAVING count >= 100
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
tsv-sort -k2,2 \
> L1.tsv

# 含有菌种id、并且属内组装在'Complete Genome', 'Chromosome'水平的物种数量大于100.

# ============================================================================



for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
echo "
    SELECT
        species_id,
        species,
        COUNT(*) AS count
    FROM ar
    WHERE 1=1
        AND genus_id = ${RANK_ID}
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    GROUP BY species_id
    HAVING count >= 100
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
tsv-sort -k2,2 \
> L2.tsv

# 属内在物种水平上组装在'Complete Genome', 'Chromosome'水平的数量大于100.

# ============================================================================


for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
echo "
    SELECT
        species_id,
        species,
        COUNT(*) AS count
    FROM ar
    WHERE 1=1
        AND genus_id = ${RANK_ID}
        AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
        AND genome_rep IN ('Full') -- fully representative
    GROUP BY species_id
    HAVING count >= 100
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
tsv-sort -k2,2 \
> L3.tsv


# 属内在物种水平组装在'Complete Genome', 'Chromosome','Scaffold'水平的物种数量大于100.
# ============================================================================

for RANK_ID in $(cat genus.list.tsv | cut -f 1); do
echo "
    SELECT
        species_id,
        species,
        COUNT(*) AS count
    FROM ar
    WHERE 1=1
        AND genus_id = ${RANK_ID}
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    GROUP BY species_id
    HAVING count >= 2
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
tsv-sort -k2,2 \
> L4.tsv

# 属内在物种水平上组装在'Complete Genome', 'Chromosome'水平的数量大于2.
# ===========================================================================
```
<!-- #endregion -->




<!-- #region -->
<!-- 不同水平的物种信息提取情况查看以及以丁香假单胞菌为例的一些特殊情况。 -->
```shell

wc -l L*.tsv
#    3 L1.tsv    
#   43 L2.tsv
#  114 L3.tsv
# 1734 L4.tsv
# 1894 total
# 一行代表一个物种，但是一个物种可以有多个组装好的基因组




for L in L1 L2 L3 L4; do
    cat ${L}.tsv |
        tsv-summarize --sum 3
done
#817
#15299
#80300
#28732
# tsv-summarize对标签分隔值文件中的字段进行聚合操作。 --sum <字段列表>[:STR] 数值之和。(仅限数值的字段）。
#数量代表的就是有多少个基因组


# 将L3中与L2中第一列相同的行（species_id）的剔除，剩下的统计还有多少个基因组。
cat L3.tsv |
    tsv-join -f L2.tsv -k 1 -e |
    tsv-summarize --sum 3
#13133
# tsv-join将输入行（"数据流"）与来自一个'过滤器'文件的行进行匹配。匹配是基于单个字段或整个行。字段可以通过字段号或字段名来指定。
# -e --exclude 排除匹配的记录。



cat L4.tsv |
    tsv-join -f L2.tsv -k 1 -e |
    tsv-join -f L3.tsv -k 1 -e |
    tsv-summarize --sum 3


# Some species are divided into several separate species
# 一些物种被分为几个独立的物种 


cd ~/data/Bacteria/summary

nwr member Pseudomonas -r species |
    grep -v " sp." |
    grep -E "syringae|genomosp"
#251701  Pseudomonas syringae group genomosp. 3  species Bacteria
#251699  Pseudomonas syringae group genomosp. 7  species Bacteria
#317659  Pseudomonas syringae pv. coryli species Bacteria
#317     Pseudomonas syringae    species Bacteria

echo "
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND species_id IN (251701,251699,317659)
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    GROUP BY species_id
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite


```
<!-- #endregion -->


### Model organisms


<!-- #region -->
<!-- 具有参考基因组信息的物种提取，筛选条件为细菌和古菌的"Proteobacteria" 变形菌门  Pseudomonadota 假单胞菌门 -->
```shell
cd ~/data/Bacteria/summary

GENUS=$(  # 细菌和古菌所有的属名
    cat genus.list.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)
echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

#提取数据库中细菌和真菌所有属的含有参考基因组的物种（15种）
# organism_name
# Caulobacter vibrioides NA1000
# Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819
# Pseudomonas aeruginosa PAO1
# Acinetobacter pittii PHEA-2
# Escherichia coli str. K-12 substr. MG1655
# Escherichia coli O157:H7 str. Sakai
# Klebsiella pneumoniae subsp. pneumoniae HS11286
# Salmonella enterica subsp. enterica serovar Typhimurium str. LT2
# Shigella flexneri 2a str. 301
# Coxiella burnetii RSA 493
# Chlamydia trachomatis D/UW-3/CX
# Staphylococcus aureus subsp. aureus NCTC 8325
# Bacillus subtilis subsp. subtilis str. 168
# Listeria monocytogenes EGD-e
# Mycobacterium tuberculosis H37Rv

cat reference.tsv |
    sed '1s/^/#/' |
    nwr append stdin -r phylum -r class |  # 加上门和纲的信息
    tsv-select -H -f 1,2,phylum,class |
    parallel --col-sep "\t" -j 1 '
        if [[ "{3}" == "Proteobacteria" || "{3}" == "Pseudomonadota" ]]; then
            printf "%s\t%s\t%s\n" {1} {2} {4}
        else
            printf "%s\t%s\t%s\n" {1} {2} {3}
        fi
    ' |
    mlr --itsv --omd cat


cp reference.tsv ~/Scripts/genomes/assembly/Bacteria.reference.tsv
# "Proteobacteria" 变形菌门  Pseudomonadota 假单胞菌门
# parallel --col-sep 把文件中的行切分为列，做为输入参数。

# if判断的三种写法，第三种支持正则：
# test expression
# [ expression ]
# [[ expression ]]

# tsv-select 读取文件或标准输入，并将选定的字段写到标准输出。字段是按照列出的顺序写的。这类似于Unix的'cut'，但有能力对字段重新排序。-H 第一行作为表头；-f 输出格式。
# mlr --itsv --omd 转换成markdown格式的表格
```
<!-- #endregion -->


## Download all assemblies
### Create assembly.tsv


<!-- #region -->
<!-- 根据上面的四个水平的物种信息和参考物种信息创建下载物种的tsv，原则同一物种的数量以最高水平的文件为标准。 -->
<!-- reference+L2+L3(具有完整参考序列)+L4 -->
<!-- 并且对提取到的物种信息进行了菌名称的缩写 -->
```shell
cd ~/data/Bacteria/summary

# 首先将具有参考基因组的物种写入raw.tsv文件中
# ============================================================================


cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level \
    > raw.tsv



# 然后将L2中的物种写入raw.tsv中，排除未知种
# =============================================================================


# L2 属内在物种水平上组装在'Complete Genome', 'Chromosome'水平的数量大于100.
SPECIES=$(
    cat L2.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)
# SPECIES 物种id集合

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter -H --regex '2:^[A-Z]' \    #第二列species物种开头是大写字母
    >> raw.tsv


# =============================================================================


# L3 属内在物种水平组装在'Complete Genome', 'Chromosome','Scaffold'（含有完整的参考序列）水平的物种数量大于100.

SPECIES=$(
    cat L3.tsv |
        tsv-join -f L2.tsv -k 1 -e |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome', 'Scaffold')
        AND genome_rep IN ('Full') -- fully representative
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter -H --regex '2:^[A-Z]' \
    >> raw.tsv


# ============================================================================

# L4 属内在物种水平上组装在'Complete Genome', 'Chromosome'水平的数量大于2
SPECIES=$(
    cat L4.tsv |
        tsv-join -f L2.tsv -k 1 -e |
        tsv-join -f L3.tsv -k 1 -e |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome') -- complete genomes
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter -H --regex '2:^[A-Z]' \
    >> raw.tsv


datamash check < raw.tsv
#38788 lines, 5 fields
# organism_name（strain_name）,species,genus,ftp_path,assembly_level

# ============================================================================




# 因为基因组的的organism_name太长，缩写
# Create abbr. （Abbreviate strain scientific names.缩写菌种学名）
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 |
    tsv-filter -H --regex '1:^[A-Z]' \
    > Bacteria.assembly.tsv  
    #abbr. ftp_path species assembly_level





# abbr_name.pl
# Usage:
#         cat <file> | perl abbr_name.pl [options]
#           Options:
#             --help              brief help message
#             --column    -c  STR Columns of strain, species, genus, default is 1,2,3.
#                                 If there's no strain, use 1,1,2.
#                                 Don't need the strain part, use 2,2,3
#                                 When there's only strain, use 1,1,1
#             --separator -s  STR separator of the line, default is "\s+"
#             --min INT           mininal length for abbreviation of species
#             --tight             no underscore between Genus and species
#             --shortsub          clean subspecies parts









datamash check < Bacteria.assembly.tsv
#38662 lines, 4 fields

# find potential duplicate strains or assemblies 
cat Bacteria.assembly.tsv |
    tsv-uniq -f 1 --repeated
# 检查缩写没有重复

cat Bacteria.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp  
# 检查有没有第二段不包含ftp的行

cat Bacteria.assembly.tsv |
    tsv-filter --or --str-in-fld 1:genomosp --str-in-fld 1:genomovar
# 检查第一段包含genomosp和genomovar


# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Bacteria.assembly.tsv
# cp Bacteria.assembly.tsv ~/Scripts/genomes/assembly

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```
<!-- #endregion -->



### rsync and check


<!-- #region -->
<!-- 根据Bacteria.assembly.tsv进行基因组下载 -->

```shell
cd ~/data/Bacteria

cat ~/Scripts/genomes/assembly/Bacteria.assembly.tsv |
    tsv-filter -v --str-in-fld 2:http

nwr assembly ~/Scripts/genomes/assembly/Bacteria.assembly.tsv \
    -o ASSEMBLY

# Remove dirs not in the list
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    tr "/" "\t" |
    cut -f 2 |
    tsv-join --exclude -k 1 -f ASSEMBLY/url.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm -fr ASSEMBLY/{}
    '

# Run
proxychains4 bash ASSEMBLY/rsync.sh

# Check md5
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

# ==> Achr_dele_GCF_013116765_2
# md5sum: md5checksums.txt: No such file or directory
# ==> Achr_dele_GCF_016127315_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Achr_dele_GCF_021432025_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Achr_deni_GCF_001514355_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Acin_bau_GCF_019264845_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Acin_johnsonii_GCF_004337595_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Acin_johnsonii_GCF_008180305_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Acin_johnsonii_GCF_015602645_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Acin_johnsonii_GCF_016027055_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Actinomy_oris_GCF_027945475_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Actinomy_oris_K20_GCF_023169925_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Actinos_pre_GCF_002354875_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Actinos_pre_pretiosum_GCF_003516205_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Actinos_pre_pretiosum_GCF_018139085_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Acu_muris_GCF_002201475_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Acu_muris_GCF_016697365_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Ad_equol_DSM_19450_GCF_000478885_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Aerom_hydrop_GCF_022631195_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Aerom_hydrop_GCF_022700835_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Aerom_hydrop_GCF_022700855_1
# md5sum: md5checksums.txt: No such file or directory
# ==> Aerom_hydrop_GCF_022759545_1
# md5sum: md5checksums.txt: No such file or directory



# Collect
bash ASSEMBLY/collect.sh

# Temporary files, possibly caused by an interrupted rsync process
find ASSEMBLY/ -type f -name ".*" > ASSEMBLY/temp.list

# -name "*"








cat ASSEMBLY/temp.list |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        if [[ -f {} ]]; then
            echo Remove {}
            rm {}
        fi
    '

# Remove ASSEMBLY/Aerom_hydrop_GCF_022700855_1/.GCF_022700855.1_ASM2270085v1_genomic.gbff.gz.MToC2b
# Remove ASSEMBLY/Aerom_hydrop_GCF_022759545_1/.GCF_022759545.1_ASM2275954v1_genomic.gbff.gz.pX94Ql
# Remove ASSEMBLY/Aerom_hydrop_GCF_022631195_1/.GCF_022631195.1_ASM2263119v1_protein.gpff.gz.Bcm0TO
# Remove ASSEMBLY/Aerom_hydrop_GCF_022700835_1/.GCF_022700835.1_ASM2270083v1_protein.gpff.gz.lfS04j

```

<!-- #endregion -->




