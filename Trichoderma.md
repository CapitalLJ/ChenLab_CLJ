## Build alignments across a eukaryotic taxonomy rank

Genus Trichoderma as an example.

### Strain info

#### list all ranks 


There are no noteworthy classification ranks other than species.

```shell
nwr member Trichoderma |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat

nwr lineage Trichoderma |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tFungi/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

# 查看木霉属的分类上的基本情况
```





#### Species with assemblies

Check also the family Hypocreaceae for outgroups.

```shell
mkdir -p ~/genome/Fungi_Trichoderma/summary
cd ~/genome/Fungi_Trichoderma/summary

# should have a valid name of genus
nwr member Hypocreaceae -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

# 和木霉同一科的所有属


wc -l genus.list.tsv
#18 genus.list

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

# 木霉同属的所有参考物种基因组信息


cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

# 木霉同属的所有物种信息



wc -l RS*.tsv GB*.tsv
#   8 RS1.tsv
#  37 GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     8
#GB1     113

```

###  Download all assemblies

#### Create assembly.tsv
If a refseq assembly is available, the corresponding genbank one is not downloaded

```shell
cd ~/genome/Fungi_Trichoderma/summary

# 酵母菌属的参考菌株的基因组信息
echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus IN ('Saccharomyces')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv


# 木霉同科的各属的参考菌株基因组信息（恰好只有木霉属的物种）
# RS
SPECIES=$(
    cat RS1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv



# Preference for refseq
cat raw.tsv |
    tsv-select -H -f "assembly_accession" \
    > rs.acc.tsv

# 参考菌株组装基因组编号



# 木霉同科的各属的菌株基因组信息（去除参考菌株基因组）
# GB
SPECIES=$(
    cat GB1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv






cat raw.tsv |
    tsv-uniq |
    datamash check
#115 lines, 6 fields

# Create abbr.  （创建简写名称的木霉属水平的各菌株基因组下载文件）
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    tsv-select -f 1-6 |
    perl ~/Scripts/genomes/bin/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\tbiosample\tspecies\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[6]}++; # abbr_name
        $seen{$F[6]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\t%s\n}, $F[6], $F[3], $F[4], $F[1], $F[5];
        ' |
    tsv-filter --or --str-in-fld 2:ftp --str-in-fld 2:http |
    keep-header -- tsv-sort -k4,4 -k1,1 \
    > Trichoderma.assembly.tsv

datamash check < Trichoderma.assembly.tsv
#114 lines, 5 fields

# find potential duplicate strains or assemblies
cat Trichoderma.assembly.tsv |
    tsv-uniq -f 1 --repeated
# 检查有没有重复
cat Trichoderma.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp
# 检查下载链接是否正确

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Trichoderma.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Trichoderma.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv


```

#### Count before download

strains.taxon.tsv - taxonomy info: species, genus, family, order, and class

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --rank genus

# --count: Count/
#     * One TSV file
#         * species.tsv
#     * Three Bash scripts
#         * strains.sh - strains.taxon.tsv, species, genus, family, order, and class
#         * rank.sh - count species and strains
#         * lineage.sh - count strains

# species.tsv   两列：简写名称加上物种名species



# strains.taxon.tsv
bash Count/strains.sh

# cat species.tsv |
# nwr append stdin -c 2 -r genus -r family -r order -r class \
#     > strains.taxon.tsv


# genus.lst and genus.count.tsv
bash Count/rank.sh

# 统计木霉属以及各属的species和strains数量
# cat genus.lst |
#     parallel --no-run-if-empty --linebuffer -k -j 4 '
#         n_species=$(
#             cat strains.taxon.tsv |
#                 tsv-filter --str-eq "3:{}" |
#                 tsv-select -f 3,2 |
#                 tsv-uniq |
#                 wc -l
#         )

#         n_strains=$(
#             cat strains.taxon.tsv |
#                 tsv-filter --str-eq "3:{}" |
#                 tsv-select -f 3,1 |
#                 tsv-uniq |
#                 wc -l
#         )

#         printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
#     ' |
#     tsv-sort -k1,1 |
#     (echo -e 'genus\t#species\t#strains' && cat) \
#     > genus.count.tsv


mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    mlr --itsv --omd cat


```


#### Download and check
- When rsync.sh is interrupted, run check.sh before restarting
- For projects that have finished downloading, but have renamed strains, you can run reorder.sh to avoid re-downloading
  - misplaced.tsv
  - remove.list
- The parameters of n50.sh should be determined by the distribution of the description statistics
- collect.sh generates a file of type .csv, which is intended to be opened by spreadsheet software.
- finish.sh generates the following files
  - omit.lst - no annotations
  - collect.pass.csv - passes the n50 check
  - pass.lst - passes the n50 check
  - rep.lst - representative or reference strains
  - counts.tsv

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --ass

# --ass: ASSEMBLY/
#     * One TSV file
#         * url.tsv
#     * And five Bash scripts
#         * rsync.sh
#         * check.sh
#         * n50.sh [LEN_N50] [N_CONTIG] [LEN_SUM]
#         * collect.sh
#         * finish.sh

bash ASSEMBLY/rsync.sh

# cat url.tsv |
#     tsv-join -f check.lst -k 1 -e |
#     parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
#         echo >&2
#         log_info "{3}\t{1}"
#         mkdir -p "{3}/{1}"
#         rsync -avP --no-links {2}/ {3}/{1}/ --exclude="assembly_status.txt"
#     '

bash ASSEMBLY/check.sh

# 比较check.list和url.tsv来显示下载完成度和更新check.list

# cat check.lst |
#     tsv-uniq |
#     tsv-join -f url.tsv -k 1 \
#     > tmp.list
# mv tmp.list check.lst

# cat url.tsv |
#     tsv-join -f check.lst -k 1 -e |
#     parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
#         if [[ ! -e "{3}/{1}" ]]; then
#             exit
#         fi
#         log_debug "{3}\t{1}"
#         cd "{3}/{1}"
#         md5sum --check md5checksums.txt --status
#         if [ "$?" -eq "0" ]; then
#             echo "{1}" >> ../../check.lst
#         else
#             log_warn "{1} checksum failed"
#         fi
#     '

bash ASSEMBLY/n50.sh 100000 1000 1000000

# 10000：最小n50长度
# 1000： contig数量小于1000
# 1000000: 总长度大于1000000.

# cat url.tsv |
#     tsv-join -f n50.tsv -k 1 -e |
#     parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
#         if [[ ! -e "{3}/{1}" ]]; then
#             exit
#         fi
#         log_debug "{3}\t{1}"

#         find "{3}/{1}" -type f -name "*_genomic.fna.gz" |
#             grep -v "_from_" | # exclude CDS and rna
#             xargs cat |
#             faops n50 -H -S -C stdin | # do not display header
#             (echo -e "{1}" && cat) |
#             datamash transpose
#     ' \
#     > tmp1.tsv

# # Combine new results with the old ones
# cat n50.tsv tmp1.tsv |
#     tsv-uniq |
#     keep-header -- sort \
#     > tmp2.tsv
# mv tmp2.tsv n50.tsv
# rm tmp*.tsv

# # Filter results with custom criteria
# cat n50.tsv |
#     tsv-filter -H --ge "N50:${LEN_N50}" |
#     tsv-filter -H --le "C:${N_CONTIG}" |
#     tsv-filter -H --ge "S:${LEN_SUM}" |
#     tr "\t" "," \
#     > n50.pass.csv


cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"

    # GCF是NCBI中对参考基因组的缩写      

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"



# bash ASSEMBLY/collect.sh

# cat url.tsv |
#     tsv-join -f check.lst -k 1 |
#     parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
#         log_debug "{3}\t{1}"
#         find "{3}/{1}" -type f -name "*_assembly_report.txt" |
#             xargs cat |
#             perl -nl -e '\''
#                 BEGIN { our %stat = (); }

#                 m{^#\s+} or next;
#                 s/^#\s+//;
#                 @O = split /\:\s*/;
#                 scalar @O == 2 or next;
#                 $O[0] =~ s/\s*$//g;
#                 $O[0] =~ s/\W/_/g;
#                 $O[1] =~ /([\w =.-]+)/ or next;
#                 $stat{$O[0]} = $1;

#                 END {
#                     my @c;
#                     for my $key ( qw( Organism_name Taxid Assembly_name Infraspecific_name BioSample BioProject Submitter Date Assembly_type Release_type Assembly_level Genome_representation WGS_project Assembly_method Genome_coverage Sequencing_technology RefSeq_category RefSeq_assembly_accession GenBank_assembly_accession ) ) {
#                         if (exists $stat{$key}) {
#                             push @c, $stat{$key};
#                         }
#                         else {
#                             push @c, q();
#                         }
#                     }
#                     print join(q(,), q({1}), @c);
#                 }
#             '\'' \
#             >> collect.csv
#     '
# 收集已经下载的基因组的详细信息 114行，113个基因组


bash ASSEMBLY/finish.sh

# tsv-join \
#     collect.csv \
#     --delimiter "," -H --key-fields 1 \
#     --filter-file n50.pass.csv \
#     > collect.pass.csv

# cat "collect.pass.csv" |
#     sed '1d' |
#     tsv-select -d, -f 1 \
#     > pass.lst

# log_info "Strains without protein annotations"
# cat url.tsv |
#     parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
#         if ! compgen -G "{3}/{1}/*_protein.faa.gz" > /dev/null; then
#             echo {1}
#         fi
#         if ! compgen -G "{3}/{1}/*_cds_from_genomic.fna.gz" > /dev/null; then
#             echo {1}
#         fi
#     ' |
#     tsv-uniq \
#     > omit.lst   #查看是否下载了蛋白序列信息或者CDS序列信息

# log_info "Representative or reference strains"
# cat collect.pass.csv |
#     tsv-filter -H -d, --not-empty "RefSeq_category" |
#     tsv-select -H -d, -f name |
#     sed '1d' \
#     > rep.lst   # 下载的基因组含有注释基因组

# log_info "Counts of lines"
# printf "#item\tcount\n" \
#     > counts.tsv

# for FILE in \
#     url.tsv check.lst collect.csv \
#     n50.tsv n50.pass.csv \
#     collect.pass.csv pass.lst \
#     omit.lst rep.lst \
#     ; do
#     cat ${FILE} |
#         wc -l |
#         FILE=${FILE} perl -nl -MNumber::Format -e '
#             printf qq($ENV{FILE}\t%s\n), Number::Format::format_number($_, 0,);
#             ' \
#         >> counts.tsv
# done

# 统计下载的文件行数：
# collect.pass.csv ： 下载的基因组通过了n50的检验；pass.lst 取第一列名称。
# omit.lst ：没有蛋白序列信息，文件只有一列，基因组名 81个基因组
# rep.list : 含有参考序列的基因组 31个基因组






cp ASSEMBLY/collect.pass.csv summary/


cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'


```

#### Rsync to hpcc

```shell

rsync -avP \
    ~/data/Trichoderma/ \
    wangq@202.119.37.251:data/Trichoderma

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Trichoderma/ \
    wangq@58.213.64.36:data/Trichoderma

# rsync -avP wangq@202.119.37.251:data/Trichoderma/ ~/data/Trichoderma

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Trichoderma/ ~/data/Trichoderma


```

#### BioSample

```shell

cd ~/data/Trichoderma

ulimit -n `ulimit -Hn`

nwr template ~/genome/Fungi_Trichoderma/assembly/Trichoderma.assembly.tsv \
    --bs

# * --bs: BioSample/
#     * One TSV file
#         * sample.tsv
#     * And two Bash scripts
#         * download.sh
#         * collect.sh [N_ATTR]

head BioSample/sample.tsv
# SAMD00028324   T_atrov_JCM_9410_GCA_001599035_1        Trichoderma_atroviride
# SAMD00028335   T_koningii_JCM_1883_GCA_001950475_1     Trichoderma_koningii
# SAMD00235762   T_asperellum_IC_1_GCA_013423425_1       Trichoderma_asperellum


bash BioSample/download.sh

# cat sample.tsv |
#     parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
#         mkdir -p "{3}"
#         if [ ! -s "{3}/{1}.txt" ]; then
#             log_debug "{1}\t{3}\t{2}"
#             curl -fsSL "https://www.ncbi.nlm.nih.gov/biosample/?term={1}&report=full&format=text" -o "{3}/{1}.txt"
#         fi
#     '


# Ignore rare attributes
bash BioSample/collect.sh 10
# 下载信息 在下载的biosample中的txt文件中提取.


datamash check < BioSample/biosample.tsv
#111 lines, 37 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/



```

#### MinHash
Estimate nucleotide divergences among strains.

Abnormal strains

If the maximum value of ANI between strains within a species is greater than 0.05, the median and maximum value will be reported. Strains that cannot be linked by the median ANI, e.g., have no similar strains in the species, will be considered as abnormal strains.
It may consist of two scenarios:
Wrong species identification
Poor assembly quality
Non-redundant strains

If the ANI value between two strains within a species is less than 0.005, the two strains are considered to be redundant.
Need these files: representative.lst and omit.lst
MinHash tree

A rough tree is generated by k-mean clustering.
These abnormal strains should be manually checked to determine whether to include them in the subsequent steps.

```shell

cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --mh \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005 \
    --height 0.4

# --mh: MinHash/
#     * One TSV file
#         * species.tsv
#     * And five Bash scripts
#         * compute.sh
#         * species.sh
#         * abnormal.sh
#         * nr.sh
#         * dist.sh

head MinHash/species.tsv

# 基因组名加物种名
# C_pro_CCMJ2080_GCA_004303015_1  Cladobotryum_protrusum
# E_web_EWB_GCA_003055145_1       Escovopsis_weberi
# E_web_GCA_001278495_1   Escovopsis_weberi
# H_perniciosus_HP10_GCA_008477525_1      Hypomyces_perniciosus





# Compute assembly sketches
bash MinHash/compute.sh


# cat species.tsv |
# tsv-join -f ../ASSEMBLY/pass.lst -k 1 |   #n50质量合格的基因组列表
# parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
#         if [[ -e "{2}/{1}.msh" ]]; then
#             exit
#         fi
#         log_info "{2}\t{1}"
#         mkdir -p "{2}"

#         find ../ASSEMBLY/{2}/{1} -name "*_genomic.fna.gz" |
            # grep -v "_from_" |  # _rna_from_genomic.fna.gz（cds_from）
#             xargs gzip -dcf |
#             mash sketch -k 21 -s 100000 -p 2 - -I "{1}" -o "{2}/{1}"
#     '
# 创建msh文件





# Distances within species
bash MinHash/species.sh

# cat species.tsv |
# tsv-join -f ../ASSEMBLY/pass.lst -k 1 |
# tsv-select -f 2 |
#     tsv-uniq \
#     > species.lst

# cat species.lst |
# while read SPECIES; do
#     log_debug "${SPECIES}"
#     mkdir -p "${SPECIES}"

#     cat species.tsv |
# tsv-join -f ../ASSEMBLY/pass.lst -k 1 |
# tsv-filter --str-eq "2:${SPECIES}" |
#         tsv-select -f 1 \
#         > "${SPECIES}/assembly.lst"      同一个物种的不同菌株基因组组装

# #    echo >&2 "Number of assemblies >= 2"
#     N_ASM=$(
#         cat "${SPECIES}/assembly.lst" | wc -l
#     )
#     if [[ $N_ASM -lt 2 ]]; then
#         continue
#     fi
    # 同一个物种含有多个基因组组装，进行下一步


#     echo >&2 "    mash distances"
#     if [[ ! -s "${SPECIES}/mash.dist.tsv" ]]; then
#         mash triangle -E -p 8 -l <(
#             cat "${SPECIES}/assembly.lst" |
#                 parallel --no-run-if-empty --linebuffer -k -j 1 "
#                     if [[ -e ${SPECIES}/{}.msh ]]; then
#                         echo ${SPECIES}/{}.msh
#                     fi
#                 "
#             ) \
#             > "${SPECIES}/mash.dist.tsv"
#     fi

    # 计算同一个物种含有多个基因组之间的遗传距离
# done




# Abnormal strains 异常菌株
bash MinHash/abnormal.sh

# ANI_VALUE_THRESHOLD=0.05

# log_info Abnormal strains

# cat species.lst |
# while read SPECIES; do
# #    log_debug "${SPECIES}"

#     # Number of assemblies >= 2
#     if [[ ! -s "${SPECIES}/mash.dist.tsv" ]]; then
#         continue
#     fi

#     D_MAX=$(
#         cat "${SPECIES}/mash.dist.tsv" |
#             tsv-summarize --max 3
#     )
#     if (( $(echo "$D_MAX < $ANI_VALUE_THRESHOLD" | bc -l) )); then
#         continue                     
#     fi              #同物种内不同菌株之间的最大遗传距离小于0.05，继续往下，否则跳此循环，继续进行下一个物种


#     # "Link assemblies with median ANI"
#     D_MEDIAN=$(
#         cat "${SPECIES}/mash.dist.tsv" |
#             tsv-filter --lt "3:$ANI_VALUE_THRESHOLD" |
#             tsv-summarize --median 3
#     )
#     cat "${SPECIES}/mash.dist.tsv" |
#         tsv-filter --ff-str-ne 1:2 --le "3:$D_MEDIAN" |
#         perl -nla -F"\t" -MGraph::Undirected -e '
#             BEGIN {
#                 our $g = Graph::Undirected->new;
#             }

#             $g->add_edge($F[0], $F[1]);

#             END {
#                 for my $cc ( $g->connected_components ) {
#                     print join qq{\n}, sort @{$cc};
#                 }
#             }
#         ' \
#         > "${SPECIES}/median.cc.lst"   遗传距离小于中位数的基因组之间通过无向图找出连通的点。

#     log_info "${SPECIES}\t${D_MEDIAN}\t${D_MAX}"
#     cat ${SPECIES}/assembly.lst |
#         grep -v -Fw -f "${SPECIES}/median.cc.lst"
# done |
#     tee abnormal.lst




cat MinHash/abnormal.lst
#T_har_CGMCC_20739_GCA_019097725_1
#T_har_Tr1_GCA_002894145_1
#T_har_ZL_811_GCA_021186515_1




# Non-redundant strains within species
bash MinHash/nr.sh


# ANI_VALUE_THRESHOLD=0.005
# log_info Non-redundant strains

# cat species.lst |
# while read SPECIES; do
#     log_debug "${SPECIES}"

#     # Number of assemblies >= 2
#     if [[ ! -s "${SPECIES}/mash.dist.tsv" ]]; then
#         continue
#     fi

#     echo >&2 "    List NR"
#     cat "${SPECIES}/mash.dist.tsv" |
#         tsv-filter --ff-str-ne 1:2 --le "3:${ANI_VALUE_THRESHOLD}" \
#         > "${SPECIES}/redundant.dist.tsv"   #遗传距离小于0.005，视为冗余基因组

#     echo >&2 "    Connected components"
#     cat "${SPECIES}/redundant.dist.tsv" |
#         perl -nla -F"\t" -MGraph::Undirected -e '
#             BEGIN {
#                 our $g = Graph::Undirected->new;
#             }

#             $g->add_edge($F[0], $F[1]);

#             END {
#                 for my $cc ( $g->connected_components ) {
#                     print join qq{\t}, sort @{$cc};
#                 }
#             }
#         ' \
#         > "${SPECIES}/connected_components.tsv"  冗余基因组之间通过无向图找出连通的点。

#     echo >&2 "    Scores based on rep.lst, omit.lst, and assembly_level"
#     # score.tsv

#     cat "${SPECIES}/connected_components.tsv" |
#         perl -nla -MPath::Tiny -F"\t" -e '
#             BEGIN {
#                 our %rep = map { ($_, 1) } path( q(../ASSEMBLY/rep.lst) )->lines({chomp => 1});
#                 our %omit = map { ($_, 1) } path( q(../ASSEMBLY/omit.lst) )->lines({chomp => 1});
#             }

#             # Representative strains are preferred
#             if ( grep { $rep{$_} } @F ) {
#                 @F = grep { ! $rep{$_} } @F
#             }
#             else {
#                 shift @F;
#             }
#             printf qq(%s\n), $_ for @F;
#             ' \
#         > "${SPECIES}/redundant.lst"

#     cat "${SPECIES}/assembly.lst" |
#         tsv-join --exclude -f "${SPECIES}/redundant.lst" \
#         > "${SPECIES}/NR.lst"  非冗余列表

# done




# Distances between all selected sketches, then hierarchical clustering
bash MinHash/dist.sh

# mash triangle -E -p 8 -l <(
#     cat species.tsv |
# tsv-join -f ../ASSEMBLY/pass.lst -k 1 |
# parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
#             if [[ -e "{2}/{1}.msh" ]]; then
#                 echo "{2}/{1}.msh"
#             fi
#         '
#     ) \
#     > mash.dist.tsv     
    
#  # 所有物种的（n 50合格）的遗传距离，下三角矩阵

# log_info Fill distance matrix with lower triangle
# tsv-select -f 1-3 mash.dist.tsv |
#     (tsv-select -f 2,1,3 mash.dist.tsv && cat) |
#     (
#         cut -f 1 mash.dist.tsv |
#             tsv-uniq |
#             parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
#         cat
#     ) \
#     > mash.dist_full.tsv
# #  全矩阵


# log_info "Clusting via R hclust(), and grouping by cutree(h=0.4)"
# cat mash.dist_full.tsv |
#     Rscript -e '
#         library(readr);
#         library(tidyr);
#         library(ape);

#         pair_dist <- read_tsv(file("stdin"), col_names=F);
#         tmp <- pair_dist %>%
#             pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
#         tmp <- as.matrix(tmp)
#         mat <- tmp[,-1]
#         rownames(mat) <- tmp[,1]

#         dist_mat <- as.dist(mat)
#         clusters <- hclust(dist_mat, method = "ward.D2")
#         tree <- as.phylo(clusters)
#         write.tree(phy=tree, file="tree.nwk")

#         group <- cutree(clusters, h=0.4)
#         groups <- as.data.frame(group)
#         groups$ids <- rownames(groups)
#         rownames(groups) <- NULL
#         groups <- groups[order(groups$group), ]
#         write_tsv(groups, "groups.tsv")
#     '


```

#### Condense branches in the minhash tree

```shell
mkdir -p ~/data/Trichoderma/tree
cd ~/data/Trichoderma/tree

nw_reroot ../MinHash/tree.nwk Sa_cer_S288C |    #将Sa_cer_S288C节点作为新的根节点进行重新定位
    nw_order -c n - \   # nw_order命令用于重新排列树的节点顺序。参数-c n指定按照节点名称的字母顺序进行排序，而-表示从标准输入读取树的输入。
    > minhash.reroot.newick

ARRAY=(
#    'order::5'
#    'family::4'
#    'genus::3'
    'species::2'
)

rm minhash.condensed.map
CUR_TREE=minhash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/genome/Scripts/condense_tree.sh ${CUR_TREE} ../Count/strains.taxon.tsv 1 ${GROUP_COL}
#  Bash 脚本，用于压缩（condense）一个 Newick 格式的进化树



    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.species.newick |
    rsvg-convert -o Trichoderma.minhash.png




```

#### Count valid species and strains for protein families

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus



# strains.taxon.tsv
bash Count/strains.sh

# genus.lst genus.count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    mlr --itsv --omd cat

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv


# 通过基因组装配水平的检查、参考序列是否完、同一物种间不同菌株间基因组是否冗余、同一物种间不同菌株遗传距离过大等筛选后剩下的菌株基因组。


# genus	#species	#strains
# Escovopsis	1	1
# Saccharomyces	1	1
# Trichoderma	16	24
```

#### Collect proteins

```shell

cd ~/data/Trichoderma/

nwr template ~/genome/Fungi_Trichoderma/summary/Trichoderma.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst

# * --pro: Protein/
#     * One TSV file
#         * species.tsv
#     * collect.sh

# collect proteins
bash Protein/collect.sh

cat Protein/counts.tsv |
    mlr --itsv --omd cat

```
| #item | count |
| --- | --- |
| Proteins | 275,985 |
| Unique headers and annotations | 275,985 |
| Unique proteins | 275,985 |
| all.replace.fa | 275,985 |
| all.annotation.tsv | 275,986 |
| all.info.tsv | 275,986 |


#### Phylogenetics with fungi61
Find corresponding proteins by hmmsearch

```shell
cd ~/data/Trichoderma

# The fungi61 HMM set
nwr kb fungi61 -o HMM

E_VALUE=1e-20

# Find all genes
for marker in $(cat HMM/fungi61.lst); do
    echo >&2 "==> marker [${marker}]"

    mkdir -p Protein/${marker}

    cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/hmm/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/${marker}/replace.tsv

    echo >&2
done

```

#### Align and concat marker genes to create species tree

```shell
cd ~/data/Trichoderma

cat HMM/fungi61.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#25      25      56

cat HMM/fungi61.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:20 --le 2:30 |
    cut -f 1 \
    > Protein/fungi61.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat HMM/fungi61.lst |
    grep -v -Fx -f Protein/fungi61.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        cat Protein/{}/replace.tsv \
            > Protein/{}/{}.replace.tsv

        faops some Protein/all.uniq.fa.gz <(
            cat Protein/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > Protein/{}/{}.pro.fa
    '

# Align each markers with muscle
cat HMM/fungi61.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"
        if [ ! -s Protein/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Protein/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in Protein/{}/{}.pro.fa -out Protein/{}/{}.aln.fa
    '

for marker in $(cat HMM/fungi61.lst); do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat Protein/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s Protein/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > Protein/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat HMM/fungi61.lst); do
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 Protein/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > Protein/fungi61.aln.fas

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f MinHash/abnormal.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
    cut -f 1 |
    fasops concat Protein/fungi61.aln.fas stdin -o Protein/fungi61.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Protein/fungi61.aln.fa -out Protein/fungi61.trim.fa -automated1

faops size Protein/fungi61.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#28706
#20432

# To make it faster
FastTree -fastest -noml Protein/fungi61.trim.fa > Protein/fungi61.trim.newick

nw_reroot Protein/fungi61.trim.newick Sa_cer_S288C |
    nw_order -c n - \
    > Protein/fungi61.reroot.newick

```

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$(
        cat groups.tsv |
            tsv-filter --str-eq 2:${TARGET_NAME} |
            tsv-select -f 1
    )
    cat groups.tsv |
        tsv-filter --str-eq 1:${SERIAL} |
        tsv-select -f 2 \
        > ${GROUP_NAME}

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> group_target.tsv

done



ARRAY=(
    'C_E_H::H_ros_CCMJ2808_GCA_011799845_1'
    'T_har_vire::T_vire_Gv29_8_GCA_020647635_1'
    'T_asperellum_atrov::T_atrov_IMI_206040_GCA_000171015_2'
    'T_lon_ree::T_ree_QM6a_GCA_000167675_2'
)


SERIAL=100
for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$((SERIAL + 1))
    GROUP_NAME_2=$(echo $GROUP_NAME | tr "_" " ")

    if [ "$GROUP_NAME" = "Trichoderma" ]; then
        cat ../ASSEMBLY/collect.pass.csv |
            tsv-filter -H -d"," --not-blank RefSeq_category |
            sed '1d' |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
        echo "C_pro_GCA_004303015_1" >> ${GROUP_NAME}
        echo "E_web_GCA_003055145_1" >> ${GROUP_NAME}
        echo "H_per_GCA_008477525_1" >> ${GROUP_NAME}
        echo "H_ros_GCA_011799845_1" >> ${GROUP_NAME}
    else
        cat ../ASSEMBLY/collect.pass.csv |
            cut -d"," -f 1,2 |
            grep "${GROUP_NAME_2}" |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
    fi

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> group_target.tsv

done

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ../../Scripts/condense_tree.sh ${CUR_TREE} ../Count/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done






















































