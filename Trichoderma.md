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


# 木霉同科的各属的参考菌株基因组信息
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


















































